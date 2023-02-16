# ADR-0003: Subselect Magma

## Status

Accepted

## Context

This ADR presents a more detailed explanation of [issue #937](https://github.com/mountetna/monoetna/issues/937), and presents our approach for addressing the performance issue raised there.

Magma `/query` and `/retrieve` API endpoints currently use SQL `JOIN` statements to collect data across multiple tables, and then chooses the columns to `SELECT` from, based upon the user's request. However, we are running into performance issues with this query architecture, as more data becomes stored in Magma. The `JOIN` architecture generates a single, virtual table in SQL with all of the requested data (each row containing some unique subset of the final query answer), and then uses Ruby code in the various `Predicate` classes to execute actions like `first`, `count`, and aggregation across the relevant rows.

For some reason, trying to stand up record instances for all these rows causes Ruby (`sequel` gem?) to sharply increase in memory consumption. In production, we can track the Magma container resource usage in Portainer, and we'll see the container killed due to an out-of-memory error (the `magma_app` service is capped at 8 GB of memory, currently). We also see a huge amount of network IO coming from the Postgres container.

Not much information can be found online regarding `sequel` memory issues. There are some threads in the `sequel-talk` Google Group that indicate it might be issues with the `DateTime` properties on our record models, i.e. `created_at` and `updated_at`. Without clear guidance on how we can tweak `sequel` performance, we decided to re-architect Magma querying instead.

This performance issue has become a blocker for implementing a Materialize workflow using `project` as the root model. A Materialize workflow would be able to automate the download of all data from a project by walking the model tree and serializing Magma records to JSON, but must first request the `project` record.

### Queries with Performance Issues

This performance limitation only impacts queries to the `/retrieve` and `/query` APIs where many one-to-many relationships exist, with many child records in each child model. Using the Magma test fixtures as examples, we would not see these performance issues when the requested predicates all exist on the `StartPredicate` model, because all of the data requested are from attributes on the main query model or a parent model and no `JOINS` are required. Here are some example queries with no performance issues. Equivalent `/retrieve` request payloads would also be okay:

- [ "monster", "::all", "stats" ]
- [ "victim", "::all", [ "birthday", "weapon", [ "monster", "name" ] ] ]

However, when we request data that traverses one-to-many relationships, then we run into this performance issue. Equivalent `/retrieve` request payloads that request all attributes / identifiers of a model with many one-to-many relationships would also reveal performance issues.

- [ "monster", "::all", [ [ "victim", "::all", [ "sidekick", "::all", "name" ] ] ] ]

### Example IPI Problem Query

This section will specifically explore the issue with IPI, though the same experience applies to all other projects. The IPI `project` model has several collection-type child models (i.e. `document`, `sc_rna_seq_pool`, `cytof_pool`, etc.). As those models get more records input, the `JOIN` architecture causes a query or retrieve action against the `project` model to explode exponentially, when requesting all attribute values (i.e. all `sc_rna_seq_pool` `tube_names`).

When running `EXPLAIN ANALYZE` in Postgres, we can see that a single query against the IPI project with all attributes results in about 1.1 million rows (as of January 2023), due to all the `JOIN` statements. This query attempts to fetch all of the child model identifiers -- effectively generating a giant table with all of the `document`, `sc_rna_seq_pool`, `rna_plate`, `cytof_pool`, etc., records that are attached to the IPI project.

Tracking the memory usage when testing locally, we see about 22-30GB of memory being consumed for the above query. Local testing also showed approximately 106 seconds to return data from a simple IPI project retrieve-all-attributes request (with local replication of production data).

The above test used a `/query` request like below, though an equivalent `/retrieve` request payload that requests all attributes of the `project` model would also reveal performance issues.

- [ "project", "::all", [ "cytof_pool", "document", "sc_rna_seq_pool", "rna_plate", "experiment" ] ]

## Decision

We decided to replace the cross-table `JOINS` with subselects. Let's examine the difference using an example query that fetches all victim names for every monster record. Note that there is a one-to-many relationship between `monster` and `victim` (one monster can have 0 or more victims):

```
[ "monster", "::all", [ [ "victim", "::all", "name" ] ] ]
```

Using the existing `JOIN` architecture, this gets converted into SQL with a `JOIN` statement between the `monster` table and the `victim` table:

```
SELECT m.monster_name as monster_name, v.victim_name as victim_name
FROM labors.monster as m,
INNER JOIN labors.victim as v
ON v.monster_id = m.id
```

Which returns a SQL result like:

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

You can imagine that if there were many one-to-many relationships stemming off of the `monster` model, this SQL result would explode and generate many rows, since each row must be a unique combination of all records from all child tables off of `monster`. For example, given three child tables with `n`, `m`, and `o` number of records, all with a single parent monster record, a query against that singular monster record would result in a table with `n * m * o` rows in it.

Using subselects instead of `JOINS`, our query gets converted into something like the following (caveat: this example, generated SQL may not exactly match the implementation):

```
SELECT m.monster_name, coalesce((
    SELECT json_agg(
        json_build_array(
            inner_m.monster_name,
            coalesce((
                SELECT json_agg(
                    v.victim_name,
                    v.victim_name
                )
                FROM labors.victim as v
                WHERE v.monster_id = inner_m.id
            ), CAST('[]' AS json))
        )
    ) FROM labors.monster as inner_m
    WHERE inner_m.id = m.id
), CAST('[]' AS json)) as victim_name
FROM labors.monster as m
```

Which generates a SQL result like this, for an `::all` verb:

| monster_name | victim_name                                                                            |
| ------------ | -------------------------------------------------------------------------------------- |
| monster1     | [ [ monster1, [ [ victim1, victim1 ], [ victim2, victim2 ], [ victim3, victim3 ] ] ] ] |
| monster2     | [ [ monster2, [ [ victim4, victim4 ], [ victim5, victim5 ] ] ] ]                       |

And this SQL result for a `::first` verb:

| monster_name | victim_name                   |
| ------------ | ----------------------------- |
| monster1     | [ victim1, victim2, victim3 ] |

One side-effect of this change is that using subselects moves the `first`, `count`, and aggregation actions into SQL, instead of doing them in Ruby, which should be more performant overall. All the answer extraction methods were changed accordingly.

Note the presence of the identifier and nested tuples for `::all` verbs. This is especially important for nested results (i.e. traversing multiple one-to-many relationships. One example is `monster` -> `victim` -> `sidekick` in Magma testing parlance), as without the identifier and tuple format, you would lose the intermediate model's identifiers in the nested results.

### Discarded SQL Aggregation Alternative

We also attempted to aggregate the `JOIN` results directly in SQL, using `json_agg` on the `JOIN` table. Something like:

```
SELECT m.monster_name as monster_name,
    json_agg(v.victim_name) as victim_name
FROM labors.monster as m,
INNER JOIN labors.victim as v
ON v.monster_id = m.id
GROUP BY m.monster_name
ORDER BY m.monster_name
```

To try and get a result like:

| monster_name | victim_name                   |
| ------------ | ----------------------------- |
| monster1     | [ victim1, victim2, victim3 ] |

However, in practice this approach was measurably less performant than the pure `JOIN` architecture. We believe this was because:

- All non-aggregated columns (i.e. the ones from the main query table / model) have to be included in the `GROUP BY` clause.
- With a large result set, the `GROUP BY` and `ORDER BY` clauses caused all the data to be written to disk, and then grouped / ordered there instead of in memory.
- The above operations were done in addition to the base `JOINS`, so we were just adding work on top of an inefficient query.

### Example IPI Problem Query

While testing locally, the same IPI project query returned in 34ms when using subselects -- a significant improvement over the `JOIN` architecture (106s).

### Subselect Classes

This ADR only affects Predicates that request data from a model different than the model of the `StartPredicate` (and typically are wrapped in a `TablePredicate`), so we introduce a new set of `Subselect` classes that are subclasses of the original `Predicate` classes. Both sets of predicates are needed to support different queries, depending on which models are accessed and their relationship to the main query model.

### Answer Classes

We implemented two main types of `Answer` classes, to accomodate `::all` and `::first` in the `ModelPredicate`.

- `AnswerTuple` is for `::all`, where we have to return a tuple, as `(<identifier>, <data>)`.
- `Answer` is for `::first`, where we just return the `<data>`.

Note that other, non-Model Predicates (i.e. all `ColumnPredicates`) that return data from the database also use the `Answer` class to encapsulate their data. This is especially important for array-like answers, such as `FileCollection`, to differentiate the result from nested `AnswerTuples`.

We also provide a `NilAnswer` class for consistency and answer formatting, to wrap `nil` from the database.

## Consequences

This change requires updates to how data is unpacked and returned from a `Magma::Question` -- now the `Predicates` have to traverse through nested arrays in the SQL response, instead of constructing responses in Ruby code. So this change touches almost every aspect of the `Question` / `/query` / `/retrieve` workflows, and it may break subtle behaviors in those actions. We have already found several edge-case bugs that were not caught by the current test suite, and we should continue to keep an eye out for more.

### Visible Improvements

- Loading the main Timur page for a project is noticeably faster. This page relies on a `/retrieve` request to get a project record and all attributes, and had significant loading time. Now, the page loads very quickly.
- We can now implement a Materialize workflow using the `project` model as a root, to fetch all data for a given project. Some IPI consortium members are exploring this tool, built into the `etna` gem (version `0.1.47` and newer).

### Remaining Inefficiencies

When debugging [GH issue #1199](https://github.com/mountetna/monoetna/issues/1199) (or perhaps a similar issue), I recall seeing that some subquery `JOINS` were still left in the main query -- I don't recall the details, but it seemed like they might not have been correctly ignored during construction of the main query. Subqueries are queries using `::any` or `::every` verbs, resulting in a SQL `HAVING` clause plus a `JOIN`. This doesn't seem to affect the subselects, but it may add some inefficiencies, especially if there wind up being a lot of subqueries. Probably something that requires further investigation and testing.

### Answer Classes

Since formatting of Answers is now encapsulated in their own classes, we can imagine that formatting for different answer environments can be more flexible. This can be seen a little bit in the `to_json` and `to_s` methods, which handle formatting the Answers for a JSON payload vs a TSV payload.

### Removal of Matrix Caching

One feature removed from Magma is the caching of `MatrixPredicate` results. The implementation cached matrix data (an array of 58k values) based on the matrix attribute's model's identifier, in a key-value store in memory. This took advantage of the fact that with the `JOIN` architecture, the matrix attribute's model columns were all available to the main `SELECT` statement, and the `MatrixPredicate` Ruby code could always access the matrix model's identifier column.

With subselects, the `MatrixPredicate` loses the ability to `SELECT` columns that are not explicitly specified by the query. This means that we do not have access to the matrix model's identifier value, and we cannot use that as a unique key in a key-value store. With our nested results, we also do not have any other unique way to identify the matrix results (any given model identifier could point to multiple matrix results, given the intervening one-to-many relationships that could exist between the `StartPredicate` model and the matrix model). Thus, we decided to remove the caching from the `MatrixPredicate` and `MatrixAttribute`, and just rely on the database results each time.
