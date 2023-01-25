---
layout: default
---

# Etna Client Gem

{:.no_toc}

The `etna` gem provides command-line clients to interact with Etna services.

- TOC
  {:toc}

## Etna Gem Usage

The `etna` gem includes several command-line utilities for managing projects and data.

### Installation

To install the etna gem, you will need a Ruby installation -- currently we run `2.7`. We recommend that you use [rbenv](https://github.com/rbenv/rbenv) to manage this. Then you can simply run:

```bash
$ gem install etna
```

Immediately after installation, you will want to do two things:

1. Source the `etna.completion` script in your shell.
2. Configure the gem's service endpoints.

#### Sourcing the Completion Script

The gem provides a completion script in your root directory (`~/etna.completion`) that allows you to tab-complete all of the commands and flags. You will need to source this script into your shell. For example, if you use `bash`, you can add the following line to your `~/.bashrc` file:

```bash
source ~/etna.completion
```

And then reload your shell (`$ source ~/.bashrc`) or open a new terminal window.

When you tab, environment variables appear with double-hyphens in front (i.e. `--environment`), and argument placeholders appear with double-underscores on both sides (i.e. `__filepath__`). **Note** that the `--environment` flag is always optional, and it defaults to the last service tier you have configured.

A full command prompt may look like:

```bash
$ etna administrate models add __project_name__ --environment
```

#### Gem Services Configuration

The gem needs to know the right service hostnames to operate against. To configure this, you need to export your Janus token into the environment and run a single command.

First, set up `etna` for production:

```bash
$ export TOKEN=<token from Janus>
$ etna config set https://<polyphemus hostname>
```

For the data science team, you should next set up the staging services. You will need to be on the VPN, grab a Janus staging token, and add an `--ignore-ssl` flag after the hostname:

```bash
$ export TOKEN=<token from Janus staging>
$ etna config set https://<polyphemus staging hostname> --ignore-ssl
```

This will make the `staging` environment your default environment to run commands against.

## Project Definition and Modeling

Two main workflows are available for defining your project models and structure, under the `administrate` command in the etna gem. These workflows are:

- Add Models: Allows you to add models and defines attributes to an existing project. Contact the engineering team to have your project created in Janus and Magma, before using this command!
- Attribute Actions: A set of four different commands that let you add, rename, or update attributes on existing models.

### Add Models

This command can be invoked with the following command:

```bash
$ etna administrate models add __project_name__
```

This will create a local CSV file that is a clone of the desired project, including models and attributes. You can use this command to download a CSV version of an existing project, and use it as a template for your project.

Running this command with a local CSV file that exists will launch a watcher that detects changes to this file and reports any validation errors in the terminal, like:

```bash
Watching for changes to mvir1_models_project_tree.csv...
Input file mvir1_models_project_tree.csv is invalid:
  * Error detected on line 314: Attribute restricted of model cytof_pool has duplicate definitions!
```

You can edit the CSV using whatever editor you like, and when you save the file the status will be updated in the terminal. Once the file passes all the validations, you will see a message like the following in the terminal:

```bash
File mvir1_models_project_tree.csv is well formatted and contains 18 models to synchronize to development mvir1.
To commit, run etna administrate models add mvir1 --file mvir1_models_project_tree.csv --target-model project --execute
```

Once you stop the validation watcher (`ctrl-c`), you can run the execute command. 

**NOTE:** Make sure your project exists in Magma and Janus before executing this command.

It will ask for you to type in a random string to verify that you want to sync your CSV to the server:

```bash
$ etna administrate models add mvir1 --file mvir1_models_project_tree.csv --target-model project --execute
File mvir1_models_project_tree.csv is well formatted and contains 18 models to synchronize to development mvir1.
Would you like to execute?
To confirm, please type LmK2hqc=:
LmK2hqc=
Executing {:action_name=>"update_attribute", :model_name=>"project", :attribute_name=>"name", :type=>nil, :description=>nil, :display_name=>"Name", :format_hint=>nil, :hidden=>nil, :index=>nil, :link_model_name=>nil, :read_only=>nil, :attribute_group=>nil, :restricted=>nil, :unique=>nil, :validation=>nil}...
...

...

Success!
```

You should now check Timur's `map` view to verify that your changes were correctly applied.

### Attribute Actions

For now, these require a JSON file input, and have a separate a command to validate the structure of your actions.

A JSON template is available for your reference. It is in the GitHub repository at the URL below, and you can also copy it out of the gem via a command like:

```bash
$ etna create_template attribute_actions
A sample attribute actions JSON template has been provided in the current directory as `attribute_actions_template.json`.
```

When using the templates, make sure to remove all the comments, which are on lines starting with `//`. JSON format does not allow comments, and they are provided for initial explanation only.

| Action            | Template                                                                                                                                     | Etna Command                               |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------ |
| Attribute Actions | [attribute_actions_template.json](https://github.com/mountetna/monoetna/blob/master/etna/lib/etna/templates/attribute_actions_template.json) | `$ etna create_template attribute_actions` |
{: rules="table"}

### Allowable Names

- Project names must be `snake_case` and not start with a number or `pg_`.
- Model names must be `snake_case` and not include numbers.
- Attribute names must be `snake_case` and not start with a number.

### Types and Required Values

#### Attributes

The following attribute keys are required for each attribute in the models:

- attribute_name
- attribute_type

The following values for `attribute_type` are supported in the CSV:

- boolean
- collection
- date_time
- file
- float
- image
- integer
- link
- match
- matrix
- string
- table

#### Models

The following model keys are required for each model\*:

- identifier
- parent_model_name
- parent_link_type

The following `parent_link_type` values are supported for model definitions:

- child
- collection
- table

**NOTE:**

- For the `project` model, **only** `identifier` is required.
- When `parent_link_type` equals `table`, **no** `identifier` key is required.

### JSON Validation

The etna gem provides commands to validate the JSON structure for the attribute actions. They are:

```bash
$ etna administrate model attributes validate_actions __project_name__ __filepath__ --environment
```

This command will output the results to the command line. Valid JSON files have a message like:

```bash
$ etna administrate models attributes validate_actions mvir1 test_attribute_actions.json
Attribute Actions JSON is well-formatted and is valid for project mvir1.
```

Whereas invalid JSON files will report a list of errors, like:

```bash
$ etna administrate models attributes validate_actions mvir1 test_attribute_actions.json
Traceback (most recent call last):
	...

  attribute_actions_from_json_workflow.rb:35:in `initialize': Attributes JSON has errors: (RuntimeError)
  * Model "assay_name" does not exist in project.
  * Model "assay_name" does not exist in project.
```

## Data Management

There is a gem command that will let you update Magma records for a single model, from a CSV file.

```bash
$ etna administrate model attributes update_from_csv __project_name__ __model_name__ __filepath__ --environment
```

The CSV file format must include column headers, with each header being the attribute name to update. The first column must be the identifier attribute for the model. You can then include rows for every record you want updated.

A simple example might be:

```bash
record_name,reference_thing,version_of_something
PROJECT001,REF001,2.7
PROJECT002,REF010,3.14
PROJECT003,REF100,9
```

**NOTE:** The command sends blank entries instead of ignoring them, so if you do not want an attribute updated for a specific record, put it into a different CSV! In the following CSV, any current data in `PROJECT002.reference_thing` will be over-written and the value set to empty string:

```bash
record_name,reference_thing,version_of_something
PROJECT001,REF001,2.7
PROJECT002,,3.14
PROJECT003,REF100,9
```

The best way to make the above update without overwriting the current value of `PROJECT002.reference_thing` would be to use two CSV files, like below:

```bash
record_name,reference_thing,version_of_something
PROJECT001,REF001,2.7
PROJECT003,REF100,9
```

```bash
record_name,version_of_something
PROJECT002,3.14
```
