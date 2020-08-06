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
types (which hold, e.g., an integer or a reference to a binary file) and link
types, which define relationships between the models in Magma.

Attached to each model is a collection of records, containing a data value for
each attribute, including links to other records and files stored on metis.

### Models

The shape of the data graph (and thus the relationships defined in each of the models) has some constraint:
- the root is always the `project` model
- every other model must have a single parent, that is, the graph is a tree.
- all links must be reciprocal, that is, for each `parent` or `link` attribute
  there must be a `child`, `collection`, or `table` attribute.
- there can be no disconnected subgraphs
- Magma strongly prefers, but does not require, a unique string identifier for each model

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

### Attributes

Each attribute has, at least, a unique attribute name within its
model, and a distinct attribute type.

The Magma attribute types are:

- {:.table-header} Value types
- `identifier` - The canonical name of a record, used to query and update the record. It is nice to have a precise validation on this attribute to ensure identifier integrity. The identifier is not strictly required; in its place Magma will substitute a database-derived id. However, Magma strongly prefers the use of identifiers, and the absence of an identifier is normally in conjunction with the `table` attribute.
- `string` - Any string.
- `integer` - Any integer (currently < 2\*\*31).
- `float` - Any floating-point number (currently does not support NaN, Infinity, or -Infinity).
- `date_time` - A date and time (usually expected in ISO8601 format)
- `file` - A reference to a file on Metis. Magma makes use of a private bucket on Metis to store files, so only Magma (not Metis) may return references to or authorize removal/update of these files.
- `image` - A reference to an image, similar to `file`.
- `matrix` - A vector of floats with fixed column labels (stored as an Array validation on the attribute). The aggregate of these vectors across the records of the model form a complete matrix, which may be sliced row-wise using record identifiers or column-wise using the column labels.


- {:.table-header} Outgoing link types:
- `parent` - Each model (except the root 'project' model) must have a single parent attribute.
- `link` - A reference to another model that is not the parent


- {:.table-header} Incoming link types:
- `child` - A one-to-one relationship.
- `collection` - A one-to-many relationship.
- `table` - A one-to-many relationship, but usually with a model lacking an identifier, making it useful to attach a table of data where the relevant identifier belongs to the parent record (e.g., a patient may have a table of drug treatments).

- {:.table-header} Other types:
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

{:.table}

### Records

What does the actual data look like? For each model, Magma can store

<div style="clear: both;"></div>

### Data restriction

### Validation

Magma models may define validations, which helps ensure the integrity of data
as it enters Magma (invalid data is rejected).  Magma has two basic forms of
validation. The first is attribute validation, which adds matchers to each
attribute on a model, e.g.:

    class MyModel < Magma::Model
      attribute :att1, type: String, match: /[r]egexp/
      attribute :att2, type: String, match: [ 'list', 'of', 'options' ]
    end

These validations are hard-coded into the model and may be hard to update. An
alternative method of validation is via a Magma::Dictionary, which allows a model
to be validated using data (records) from another model.

#### Dictionaries

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


## API

  **POST _/update_** - accepts a JSON post in this format:

    {
      "project_name" : "labors",
      "revisions" : {
        "monster" : {
          "Nemean Lion" : {
            "species" : 'lion'
          }
        }
      }
    }

  **POST _/retrieve_** - accepts a JSON post in this format:

    {
      "project_name"    : "labors",
      "model_name"      : "labor",
      "record_names"    : [ "Nemean Lion" ],
      "attribute_names" : [ "name", "number", "completed" ]
    }

  In return you will get a payload like this:

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

  **POST _/query_** - accepts JSON queries in Magma's Query language. See <https://github.com/mountetna/magma/wiki/Query> for more details.

## Setup

### Installation

We may start with a basic git checkout:

`$ git clone https://github.com/mountetna/magma.git`

Magma runs as a Ruby Magma is a Rack application, which means you can run it using any Rack-compatible server (e.g. Puma or Passenger).


### Configuration

Magma has a single YAML config file, `config.yml`; DO NOT TRACK this file, as it will hold all of your secrets. It uses the Etna::Application configuration syntax, e.g.:

    ---
    :test:
      :host: https://magma.test

    :development:
      :log_file: ./log/error.log

The environment is by default `development` but may be set via the environment variable MAGMA_ENV.

Some things you may configure:

    # This is the database configuration for the Sequel ORM (Documented at https://sequel.jeremyevans.net)
    # Magma uses the postgres adapter; it may not work with other databases.
    :db:
      :database: magma
      :host: localhost
      :adapter: postgres
      :user: magma
      :password: AAAAAAAA

    # A space-separated lists of magma project directories.
    # See below for details on creating a project.
    :project_path: ./projects/labors/

    # Where Magma will attempt to log - some server errors may not be trapped here.
    :log_file: log/error.log

    # Magma uses the `fog/aws` gem to connect to S3 for file storage,
    # and Carrierwave to manage uploads
    :storage:
      :provider: fog/aws
      :directory: 'my-magma-bucket'
      :expiration: 1440
      :credentials:
        :provider: 'AWS'
        :aws_access_key_id: 'AKIAETCETERA'
        :aws_secret_access_key: 'SoMeSecrEtK/Ey'
        :region: 'us-area-52'

    # The algorithm used by the authentication service (Janus) to sign tokens
    # and the public key to validate them.
    :token_algo: RS256
    :rsa_public: |
      -----BEGIN PUBLIC KEY-----
      KeYGoEsHeRE==
      -----END PUBLIC KEY-----

### Attributes

Attributes describe the data elements of the model, their interactions with other models, and so on.

  **_attribute_** - a generic attribute, representing any valid column name in the database

  **_identifier_** - a special attribute, a unique string identifier for records of this model type. Links between records are specified using the _identifier_ attribute. If no identifier is set the identifier defaults to the database _id_.

  <u>link types:</u>

  **_parent_** - an ancestor in a hierarchy, stored via a foriegn_key column in the model

  **_child_** - a child in a hierarchy, stored in the foreign model
  link - same as 'parent', just intended to represent other sorts of relations for clarity

  **_collection_** - an array of links to children

  <u>tables:</u>

  **_table_** - similar to a collection, except while a collection references a wholly separate model, a table is intended to essentially be a vector of values that loads with this model.

  <u>file types</u>:

  **_file_** - a generic binary document, stored on S3

  **_image_** - an image, similar to a document except it allows some thumbnailing

  <u>restriction</u>:
  **_restricted_** - a boolean attribute that determines whether this record is restricted.

### Migrations

Magma attempts to maintain a strict adherence between its models and
the database schema by suggesting migrations. These are written in the
Sequel ORM's migration language, not pure SQL, so they are fairly
straightforward to amend when Magma plans incorrectly.

To plan a new set of migrations, the first step is to amend your
models.  This also works in the case of entirely new models. Simply
sketch them out as described above, setting out the attributes each
model requires and creating links between them.

Once you've defined your models, you can execute `bin/magma plan` to
create a new migration. If you want to restrict your plan to a single
project you may do `bin/magma plan <project_name>`. Magma will output
ruby code for a migration using the Sequel ORM - you can save this in
your project's migration folder (e.g.
`project/my_project/migration/01_initial_migration.rb`).

After your migrations are in place, you can try to run them using `bin/magma
migrate`, which will attempt to run migrations that have not been run yet. If
you change your mind, you can roll backwards (depending on how reversible your
migration is) using `bin/magma migrate <migration version number>`.
