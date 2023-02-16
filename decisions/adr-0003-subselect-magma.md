# ADR-0003: Subselect Magma

## Status

Accepted

## Context

Magma `/query` and `/retrieve` API endpoints currently use SQL `JOIN` statements to collect data across multiple tables, and then defines columns to `SELECT` from. However, we are running into performance issues with this query scheme, as more data becomes stored in Magma. For example, the IPI `project` model has several collection-type child models (i.e. `document`, `sc_rna_seq_pool`, `cytof_pool`, etc.). As those models get more records input, these `JOIN` statements cause the query to explode exponentially.

When running `EXPLAIN ANALYZE` in Postgres, we can see that a single query against the IPI project with all attributes results in about 1.1 million rows (as of January 2023), due to all the `JOIN` statements. We then use Ruby code in the various `Predicate` classes, to execute actions like `first`, `count`, and aggregation across child models (foreign key relations).

For some reason, trying to stand up record instances for all these rows causes Ruby (`sequel` gem?) to sharply increase in memory consumption. Tracking the memory usage locally, we see about 22-30GB of memory being consumed. In production, we can track the Magma container resource usage in Portainer, and we'll see the container die due to an out-of-memory error (the `magma_app` service is capped at 8 GB of memory, currently). We also see a huge amount of network IO coming from the Postgres container. Local testing also showed approximately 106 seconds to return data from the IPI project query (with a simple, local replication of production data).

Not much information can be found online regarding `sequel` memory issues. There are some threads in the `sequel-talk` Google Group that indicate it might be issues with the `DateTime` properties on our record models, i.e. `created_at` and `updated_at`. Without clear guidance on how we can tweak `sequel` performance, we decided to re-architect Magma querying instead.

## Decision

We decided to replace the cross-table `JOINS` with subselects, instead. Let's examine the difference using an example query:

```
[ "monster", "::all", [ [ "victim", "::all", ":name" ] ] ]
```

Using `JOIN` statements, this was converted into SQL like:

```
SELECT m.monster_name as monster_name, v.victim_name as victim_name
FROM labors.monster as m,
INNER JOIN labors.victim as v
ON v.monster_id = m.id
```

Which generates a SQL result like:

| monster_name | victim_name |
| ------------ | ----------- |
| monster1     | victim1     |
| monster1     | victim2     |
| monster1     | victim3     |
| monster2     | victim4     |
| monster2     | victim5     |

And then in Ruby, our Predicates would handle aggregating those results into a record like:

```
{
    "monster1": {
        "victim": [ "victim1", "victim2", "victim3" ]
    },
    "monster2": {
        "victim": [ "victim4", "victim5" ]
    }
}
```

Using subselects instead of `JOINS`, our query changes into something, like:

```
SELECT m.monster_name, coalesce(SELECT json_array(json_build_array(m.monster_name, v.victim_name)) FROM labors.victim as v, CAST('[]', json)) as victim_name
FROM labors.monster as m
```

Which generates a SQL result like this, for an `::all` verb:

| monster_name | victim_name                                 |
| ------------ | ------------------------------------------- |
| monster1     | [ monster1, [ victim1, victim2, victim3 ] ] |
| monster2     | [ monster2, [ victim4, victim5 ] ]          |

And this SQL result for a `::first` verb:

| monster_name | victim_name                   |
| ------------ | ----------------------------- |
| monster1     | [ victim1, victim2, victim3 ] |

One side-effect of this change is that subselects move the `first`, `count`, and aggregation actions into SQL, instead of doing them in Ruby, which should be more performant overall. All the answer extraction methods were changed accordingly.

Note the presence of the identifier and nested tuples for `::all` verbs. This is especially important for nested results (i.e. traversing multiple one-to-many relationships. One example is monster -> victim -> sidekick in Magma testing parlance), as without the identifier and tuple format, you would lose the intermediate model's identifiers in the nested results.

While testing locally, the same IPI project query returned in 34ms -- a significant improvement over the previous impelmentation (106s).

### Magma::Answer Classes

We implemented two main types of `Answer` classes, to accomodate `::all` and `::first` in the `ModelPredicate`.

- `AnswerTuple` is for `::all`, where we have to return a tuple, as `(<record_name>, <data>)`.
- `Answer` is for `::first`, where we just return the `<data>`.

Note that other, non-Model Predicates (i.e. all `ColumnPredicates`) that return data from the database also use the `Answer` class to encapsulate their data. This is especially important for array-like answers, such as `FileCollection`, to differentiate the result from nested arrays.

We also provide a `NilAnswer` result for consistency and answer formatting, to wrap no result, or `nil` from the database.

## Consequences

This change requires updates to how data is unpacked and returned from a `Magma::Question` -- now the `Predicates` have to traverse through nested arrays in the SQL response, instead of constructing responses in Ruby code.

### Remaining Inefficiencies

When debugging GH issue #1199 (or perhaps a similar issue), I recall seeing that some subquery `JOINS` were still left in the main query -- I don't recall the details, but it seemed like they were not ignored correctly. Subqueries are queries using `::any` or `::every` verbs, resulting in a SQL `HAVING` clause plus a `JOIN`. This doesn't seem to affect the subselects, but it may add some inefficiencies, especially if there wind up being a lot of subqueries. Probably something that requires further investigation and testing.

### Magma::Answer Classes

Since formatting of Answers is now encapsulated in their own classes, we can imagine that formatting for different answer environments can be more flexible.
