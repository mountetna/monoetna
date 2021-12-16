---
layout: default
group: app
---
# Magma
{:.no_toc}
## Data Warehouse
{:.no_toc}

* TOC
{:toc}

## Organization

The purpose of Magma is to define a data graph for each Etna project: a set of
models for each of the entities in the project dataset.

Each model has a set of attributes, which may broadly be divided into value
types, which hold, e.g., an `integer` or a reference to a binary `file`, and link
types, e.g. `parent` or `collection`, which define relationships between the
models in Magma.

Attached to each model is a collection of records, containing a data value for
each attribute, including links to other records and files stored on metis.

### Models

The shape of the data graph (and thus the relationships defined in each of the models) has some constraint:
* the root is always the `project` model
* every other model must have a single parent, i.e., the graph is a tree.
* all links must be reciprocal, i.e., for each `parent` or `link` attribute
  there must be a `child`, `collection`, or `table` attribute.
* there can be no disconnected subgraphs
* Magma strongly prefers, but does not require, a unique string identifier for each model

Here is a sketch of what the graph for the "olympics" project might look like:

```
project {
  identifier project_name
  collection event
  collection athlete
}

event {
  parent project
  identifier event_name
  table entry
}

entry {
  parent event
  link athlete
  integer placement
  integer score
}

athlete {
  parent project
  identifier name
  collection entry
}

```

Here all links are reciprocal, and every model descends from the project. However,
using the `link` attribute we may indicate other one-to-one or one-to-many
relationships, which allows the graph to be more like a directed acyclic graph
(DAG) than a tree.

Models may also specify a dictionary model and mapping to be used for more
complex validations (see below on dictionary validation).

### Attributes

Each attribute has, at least, a unique attribute name within its
model, and a distinct attribute type.

The Magma attribute types are:

**Value types**
- `identifier` - The canonical name of a record, used to query and update the record. It is nice to have a precise validation on this attribute to ensure identifier integrity. The identifier is not strictly required; in its place Magma will substitute a database-derived id. However, Magma strongly prefers the use of identifiers, and the absence of an identifier is normally in conjunction with the `table` attribute.
- `string` - Any string.
- `integer` - Any integer (currently < 2\*\*31).
- `float` - Any floating-point number (currently does not support NaN, Infinity, or -Infinity).
- `boolean` - True or false.
- `date_time` - A date and time (usually expected in ISO8601 format)
- `file` - A reference to a file on Metis. Magma makes use of a private bucket on Metis to store files, so only Magma (not Metis) may return references to or authorize removal/update of these files.
- `image` - A reference to an image, similar to `file`.
- `matrix` - A vector of floats with fixed column labels (stored as an Array validation on the attribute). The aggregate of these vectors across the records of the model form a complete matrix, which may be sliced row-wise using record identifiers or column-wise using the column labels.


**Outgoing link types**
- `parent` - Each model (except the root 'project' model) must have a single parent attribute.
- `link` - A reference to another model that is not the parent


**Incoming link types**
- `child` - A one-to-one relationship.
- `collection` - A one-to-many relationship.
- `table` - A one-to-many relationship, but usually with a model lacking an identifier, making it useful to attach a table of data where the relevant identifier belongs to the parent record (e.g., a patient may have a table of drug treatments).

**Other types**
- `match` - An obscure type useful in writing dictionary-based validations
- `restricted` - a boolean value, but with special behavior that will cause Magma to censor the record from anyone lacking the `privileged` flag on their Janus token.

In addition to its type, each attribute may set several other fields:

- `description` - A pithy description of the attribute's purpose
- `display_name` - A human-legible version of the attribute name
- `format_hint` - A description of how valid contents should be formatted
- `hidden` - a flag for Whether data should be returned for this attribute
- `index` - a flag for whether to create an index on this attribute for faster querying
- `read_only` - a flag for whether the values in this attribute can be changed once set
- `restricted` - a flag for whether this attribute should be hidden from users without the `privileged` flag in their Janus token.
- `unique` - a flag for whether the values in this attribute should be unique across records.
- `validation` - a JSON value allowing validation of input values for this attribute

### Records

A record is a set of values for each attribute in the model. The set of
records for a project form a data graph with a single project record at
the root.

<div style="clear: both;"></div>

### Validation

Magma models may define validations, which helps ensure the integrity of data
as it enters Magma (invalid data is rejected).  Magma has two basic forms of
validation. The first is attribute validation, which adds matchers to each
attribute on a model, e.g.:

    class MyModel < Magma::Model
      attribute :att1, type: String, match: /[r]egexp/
      attribute :att2, type: String, match: [ 'list', 'of', 'options' ]
    end

#### Dictionaries

An alternative method of validation is via a Magma::Dictionary, which allows a
model to be validated using data (records) from another model.

You may define a dictionary relation as follows:

    class MyModel < Magma::Model
      attribute :att1
      attribute :att2

      dictionary DictModel, att1: :dict_att1, att2: :dict_att2
    end

A `my_model` record is valid if there is a matching entry in the dictionary
model, i.e., where `my_record.att1` matches `dict_record.att1` and
`my_record.att2` matches `dict_record.att2`. Here 'match' might mean
'equality', but a dictionary may also include a 'match' attribute.

    class DictModel < Magma::Model
      attribute :dict_att1
      match :dict_att2
    end

A match attribute contains json data like `{type,value}`. This allows us to construct more complex entries:

    # match any value in this Array
    { type: 'Array', value: [ 'x', 'y', 'z' ] }
    # match within this Range
    { type: 'Range', value: [ 0, 100 ] }
    # match this Regexp
    { type: 'Regexp', value: '^something$' }
    # match an ordinary value
    { type: 'String', value: 'something' }


## Usage

### Clients

The main way to interact with Magma is via the API. In the simplest case this
can be done using curl or wget to POST if one sets the
`Authorization: Etna <your token>` header; but any HTTP client will suffice.
Visit your Janus instance to get a current token, or make use of other ways to
[authenticate with Janus]({{ site.baseurl }}{% link janus.md %}#authenticating).

To use `Etna::Client` to connect to Magma in Ruby, you may `gem install etna` and then create a new client:

```
require 'etna'

e = Etna::Client.new('https://magma.example.org', ENV['TOKEN'])

payload = e.retrieve(project_name: 'labors', model_name: 'monster', record_names: [ 'Nemean Lion', 'Lernean Hydra' ], attribute_names: "all")
```

### API

The main way to interact with Magma directly is via its API (you may also perform a great many of the same operations using the data browser [Timur]({{ site.baseurl }}{% link timur.md %})).

There are four main endpoints: `update`, `retrieve`, `query`, and
`update_model`. All of them expect a POST in JSON format with a valid Etna
authorization header (i.e., `Authorization: Etna <valid janus token>`).

#### /update

**Parameters**

- `project_name` - the project key

- `revisions` - A hash with model names as keys and a set of model revisions as values. Each model revision is a hash with record identifiers as keys and a record update as value. The record update, in turn, is a hash with attribute names as keys and data as values.

  You may include any number of models, with any number of records for each model, in your revisions. Record entries need not be complete (not every attribute needs to be set), or even consistent (a different subset of attributes can be set on each record). Contrary to the name, /update will both insert new records AND update existing records.

  Valid data values vary by attribute:

	- `identifier`, `string` - a string
	- `integer` - an integer
	- `float` - a float
	- `boolean` - true or false
	- `date_time` - A date/time in an ISO8601-formatted string (e.g., "2020-02-02", or "2020-02-02T20:20")
	- `file`, `image` - A hash in the form `{ path, original_filename }`. `original_filename` is optional and may be any string. `path` may take the values `null` (to set an empty file state), `::blank` (to set a blank file state), `::temp` (to return a Metis upload url), or a valid Metis path, i.e., `metis://<project_name>/<bucket_name>/<file_path>` (to link the file into the Magma record).
	- `matrix` - An array of floats or integers.
	- `parent`, `link`, `child` - a valid string identifier for the linked model. If there is no such record, Magma will create it with the given identifier. 
	- `collection` - An array of valid string identifiers for the linked model. If there are no such records, Magma will create them with the given identifiers.
	- `table` - An array of temporary ids. The corresponding records in the
	  linked model should be inserted at the same time using the same
	  temporary id, which takes the form `::temp<any string>`.
	- `match` - A hash of the form `{type, value}`, where type may be `"Array", "Regexp", "Range", with a corresponding value.

	In addition each of these values may be set to `null`.

**Examples**

The basic revision format looks like this:

    {
      "project_name" : "labors",
      "revisions" : {
        "monster" : {
          "Nemean Lion" : {
            "species" : 'lion'
          },
	  "Lernean Hydra" : {
	    "species" : 'hydra'
	  }
        }
      }
    }

#### /retrieve

**Parameters**

Required parameters:
- `project_name` - the project key
- `model_name` - the model you wish to retrieve (in `snake_case`), or `"all"`
- `record_names` - an array of record identifiers from the model in `model_name`, or `"all"`
- `attribute_names` - an array of attribute names, or the strings `"identifier"` or `"all"`

Optional parameters:
- `collapse_tables` - whether to return records from linked models via table attributes specified in `attribute_names`
- `format` - `"json"` or `"tsv"` - the latter will force `collapse_tables: true`
- `filter` - a string defining a filter to apply to the records - a space-separated list of terms in the form `<column_name><operator><value>`, e.g. `"age>60"`. Valid operators for `string` and `identifier` columns are `=`,`~` (which will match a regular expression value) and for numeric or `date_time` columns are `=`,`>`,`<`,`<=`, `>=`
- `page_size` - Splits retrieval into sets of `page_size` records
- `page` - Retrieves page number `page` (1-indexed) records. The retrieval of page 1 will also include a count of all records.

**Examples**

A basic request for a record looks like this:

    {
      "project_name"    : "labors",
      "model_name"      : "labor",
      "record_names"    : [ "Nemean Lion" ],
      "attribute_names" : [ "name", "number", "completed" ]
    }

The output is in "payload" format, containing a hash `{ models }` keyed by model_name, and returning for each model `{ documents, template }`. The template is a complete description of the model sufficient for import into another Magma instance. The returned documents are keyed by the record identifiers, with each record containing values for the attributes requested in `attribute_names`.

    {
      "models": {
        "labor": {
          "documents": {
            "Nemean Lion": {
              "name": "Nemean Lion",
              "number": 1,
              "completed": true
            }
          },
          "template": {
            "name": "labor",
            "attributes": {
              "name": {
                "name": "name",
                "type": "String",
                "attribute_class": "Magma::Attribute",
                "display_name": "Name",
                "shown": true
              }
              // etc. for ALL attributes, not just requested
            }
          }
        },
        "identifier": "name",
        "parent": "project"
      }
    }

A few special cases exist. Here is the "template" query, which will retrieve all of the project templates but no documents:

    { "project_name": "labors", "model_name": "all", "record_names":[], "attribute_names": "all" }


The "identifier" query will retrieve all of the project identifiers at once: 

    { "project_name": "labors", "model_name": "all", "record_names": "all", "attribute_names": "identifier" }

#### /query

The Magma Query API lets you pull data out of Magma through an expressive query interface.

**Parameters**

- `project_name` - the project key
- `query` - An array of query predicates (see below for details)
- `format` - Optional. If `tsv` is provided, a TSV is returned instead of a JSON paylod. By default, a JSON payload is provided.
  - `user_columns` - If `format=tsv`, can provide an array of string values that act as column labels. The length of this array must match the number of requested columns.
  - `expand_matrices` - If `format=tsv`, boolean value to expand matrix attributes into individual columns, with one matrix data point per column.
  - `transpose` - If `format=tsv`, boolean value to transpose the resulting TSV. This means that rows become columns, and vice versa.

A general form of the query is:

    [ *predicate_args, *verb_args, *predicate_args, *verb_args, ... ]

  * predicate - an object in magma (model, record, attribute)
  * verb - a function yielding a new predicate

A basic query might look like this:

    [ 'labor', '::all', 'name' ]

A breakdown of the terms:
`labor` - specifies the model we wish to search, yielding a model predicate
`::all` - a verb argument to the model predicate, iterating across all of the items in the model, and yielding a record predicate
`name` - a verb argument to the record predicate specifying an attribute name, yielding a value and terminating the query

While the query must eventually terminate in a value (or array of values if an
array argument is passed to a record predicate), via records we might traverse
through the graph first:

    [ 'labor', '::all', 'monster', 'victim', '::first', 'city' ]

The response:

    {
       "answer" : [ [ 'Nemean Lion', 'Nemean Lion' ], [ 'Lernean Hydra', 'Lernean Hydra' ] ]
       "format" : ['labors::labor#name', 'labors::labor#name']
    }

The format describes the returned values. If the format is an array, the format
will contain a list of items with the given format. The format is usually
written in `project_name::model_name#attribute_name` format.

A more advanced query might include a filter:

    [ 'monster', [ '::has', 'stats' ], '::all', 'name' ]

Filters may be applied to any model we traverse through:

    [ 'labor', '::all', 'prize', [ 'worth', '>', '200' ], '::first', 'name' ]

**Predicates**

There are a handful of predicate types, each of which take various arguments. 

##### Model

A Model predicate is our query starting point and specifies a set of records. Model predicates can accept an arbitrary number of filter [] arguments, followed by:

    ::first - reduce this model to a single item
    ::all - return a vector of values for this model, labeled with this model's identifiers

##### Record

A Record predicate follows after a Model predicate. The valid arguments are:

    <attribute_name> - a string specifying an attribute on this model
    ::has, <attribute_name> - a boolean test for the existence of <attribute_name> (i.e., the data is not null)
    ::identifier - an alias for the attribute_name of this Model's identifier. E.g., if a Sample has identifier attribute 'sample_name', '::identifier' will return the same value as 'sample_name'

##### Column

Column attributes usually just return their value. However, you may optionally follow them with arguments to apply a boolean test.

`string`

    ::equals, <string> - A boolean test for equality, e.g. [ 'sample_name', '::equals', 'Dumbo' ]
    ::in, [ list of strings ] - A boolean test for membership, e.g., [ 'sample_name', '::in', [ 'ant', 'bear', 'cat' ] ]
    ::matches, <string> - A boolean test for a regular expression match, e.g., [ 'sample_name', '::matches', '[GD]umbo' ]

`integer`, `date_time`

    ::<= - less than or equals
    ::< - less than
    ::>= - greater than or equals
    ::> - greater than
    ::= - equals

`boolean`

    ::true - is true
    ::false - is false

`file`, `image`

    ::url - a URL to retrieve this file resource
    ::path - the filename/path for this file resource

`matrix`

    ::slice - retrieve a subset of columns from the matrix

#### /update_model

Coming soon.

## Setup

### Installation

Start with a basic git checkout:

`$ git clone https://github.com/mountetna/magma.git`

Magma is a Rack application, which means you can run it using any Rack-compatible server (e.g. Puma or Passenger).

### Configuration

Magma has a single YAML config file, `config.yml`; DO NOT TRACK this file, as it will hold all of your secrets. It uses the Etna::Application configuration syntax. See [config.yml.template](https://raw.githubusercontent.com/mountetna/magma/master/config.yml.template) for an example configuration.

### Migrations

Magma attempts to maintain a strict adherence between its models and
the database schema by suggesting migrations. These are written in the
Sequel ORM's migration language, not pure SQL, so they are fairly
straightforward to amend when Magma plans incorrectly.

To plan a new set of migrations, the first step is to amend your models.  This
also works in the case of entirely new models. Simply sketch them out as
described above, setting out the attributes each model requires and creating
links between them.

Once you've defined your models, you can execute `bin/magma plan` to create a
new migration. If you want to restrict your plan to a single project you may do
`bin/magma plan <project_name>`. Magma will output ruby code for a migration
using the Sequel ORM - you can save this in your project's migration folder
(e.g.  `project/my_project/migration/01_initial_migration.rb`).

After your migrations are in place, you can try to run them using `bin/magma
migrate`, which will attempt to run migrations that have not been run yet. If
you change your mind, you can roll backwards (depending on how reversible your
migration is) using `bin/magma migrate <migration version number>`.
