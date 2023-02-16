# ADR-0003: Subselect Magma

## Status

Accepted

## Context

This ADR is a more detailed explanation of [issue #937](https://github.com/mountetna/monoetna/issues/937), and presents our approach for addressing the performance issue raised there.

Magma `/query` and `/retrieve` API endpoints currently use SQL `JOIN` statements to collect data across multiple tables, and then defines columns to `SELECT` from. However, we are running into performance issues with this query scheme, as more data becomes stored in Magma. For example, the IPI `project` model has several collection-type child models (i.e. `document`, `sc_rna_seq_pool`, `cytof_pool`, etc.). As those models get more records input, these `JOIN` statements cause the query to explode exponentially.

When running `EXPLAIN ANALYZE` in Postgres, we can see that a single query against the IPI project with all attributes results in about 1.1 million rows (as of January 2023), due to all the `JOIN` statements. We then use Ruby code in the various `Predicate` classes, to execute actions like `first`, `count`, and aggregation across child models (foreign key relations).

For some reason, trying to stand up record instances for all these rows causes Ruby (`sequel` gem?) to sharply increase in memory consumption. Tracking the memory usage locally, we see about 22-30GB of memory being consumed. In production, we can track the Magma container resource usage in Portainer, and we'll see the container die due to an out-of-memory error (the `magma_app` service is capped at 8 GB of memory, currently). We also see a huge amount of network IO coming from the Postgres container. Local testing also showed approximately 106 seconds to return data from a simple IPI project retrieve-all-attributes request (with a local replication of production data)

Not much information can be found online regarding `sequel` memory issues. There are some threads in the `sequel-talk` Google Group that indicate it might be issues with the `DateTime` properties on our record models, i.e. `created_at` and `updated_at`. Without clear guidance on how we can tweak `sequel` performance, we decided to re-architect Magma querying instead.

### Queries with Performance Issues

This performance hit only impacts queries to the `/retrieve` and `/query` APIs where many one-to-many relationships exist, with many child records in each child model. Using the Magma test fixtures as examples, we would not see these performance issues when the requested predicates all exist on the `StartPredicate` model. The following queries would not trigger this issue because all of the data requested are from attributes on the main query model or a parent model. Equivalent `/retrieve` request payloads would also be okay:

- [ "monster", "::all", "stats" ]
- [ "victim", "::all", [ "birthday", "weapon", [ "monster", "name" ] ] ]

However, when we request data that traverses one-to-many relationships, then we run into this performance issue. Equivalent `/retrieve` request payloads that request all attributes / identifiers of a model with many one-to-many relationships (i.e. `project`) would also reveal performance issues.

- [ "monster", "::all", [ [ "victim", "::all", [ "sidekick", "::all", "name" ] ] ] ]
- [ "project", "::all", [ "cytof_pool", "document", "sc_rna_seq_pool", "rna_plate", "experiment" ] ]

## Decision

We decided to replace the cross-table `JOINS` with subselects. Let's examine the difference using an example query:

```
[ "monster", "::all", [ [ "victim", "::all", ":name" ] ] ]
```

Using `JOIN` statements, this gets converted into SQL like:

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

And then in Ruby, our Predicates would aggregate those results into a record like:

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

Using subselects instead of `JOINS`, our query gets converted into something like:

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

Note the presence of the identifier and nested tuples for `::all` verbs. This is especially important for nested results (i.e. traversing multiple one-to-many relationships. One example is `monster` -> `victim` -> `sidekick` in Magma testing parlance), as without the identifier and tuple format, you would lose the intermediate model's identifiers in the nested results.

While testing locally, the same IPI project query returned in 34ms -- a significant improvement over the previous impelmentation (106s).

### Subselect Classes

This ADR only affects Predicates that are based on a model different than the model of the `StartPredicte`, so we introduce a set of new `Subselect` classes that are subclasses of the original `Predicate` classes. Both sets of predicates are needed to support different queries.

### Answer Classes

We implemented two main types of `Answer` classes, to accomodate `::all` and `::first` in the `ModelPredicate`.

- `AnswerTuple` is for `::all`, where we have to return a tuple, as `(<record_name>, <data>)`.
- `Answer` is for `::first`, where we just return the `<data>`.

Note that other, non-Model Predicates (i.e. all `ColumnPredicates`) that return data from the database also use the `Answer` class to encapsulate their data. This is especially important for array-like answers, such as `FileCollection`, to differentiate the result from nested arrays.

We also provide a `NilAnswer` result for consistency and answer formatting, to wrap `nil` from the database.

## Consequences

This change requires updates to how data is unpacked and returned from a `Magma::Question` -- now the `Predicates` have to traverse through nested arrays in the SQL response, instead of constructing responses in Ruby code.

### Remaining Inefficiencies

When debugging GH issue #1199 (or perhaps a similar issue), I recall seeing that some subquery `JOINS` were still left in the main query -- I don't recall the details, but it seemed like they were not ignored correctly during construction of the main query. Subqueries are queries using `::any` or `::every` verbs, resulting in a SQL `HAVING` clause plus a `JOIN`. This doesn't seem to affect the subselects, but it may add some inefficiencies, especially if there wind up being a lot of subqueries. Probably something that requires further investigation and testing.

### Answer Classes

Since formatting of Answers is now encapsulated in their own classes, we can imagine that formatting for different answer environments can be more flexible. This can be seen a little bit in the `to_json` and `to_s` methods, which handle formatting the Answers for a JSON payload vs a TSV payload.

### Removal of Matrix Caching

One feature removed from Magma is the caching of `MatrixPredicate` results. The implementation cached matrix data (an array of 58k values) based on the matrix attribute's model's identifier, in a key-value store in memory. This took advantage of the fact that with the `JOIN` architecture, the matrix attribute's model columns were all available to the main `SELECT` statement, and the `MatrixPredicate` Ruby code could always access the matrix model's identifier column.

With subselects, the `MatrixPredicate` loses the ability to `SELECT` columns that are not explicitly specified by the query. This means that we do not have access to the matrix model's identifier value, and we cannot use that as a unique key in a key-value store. With our nested results, we also do not have any other unique way to identify the matrix results (any given model identifier could point to multiple matrix results, given the intervening one-to-many relationships that could exist between the `StartPredicate` model and the matrix model). Thus, we decided to remove the caching from the `MatrixPredicate` and `MatrixAttribute`, and just rely on the database results each time.
