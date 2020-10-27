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

To install the etna gem, you will need a Ruby installation -- currently we run `2.5.7`. We recommend that you use [rbenv](https://github.com/rbenv/rbenv) to manage this. Then you can simply run:

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
$ etna administrate project validate __filepath__ --environment
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

Three workflows are available for defining your project models and structure, under the `administrate` command in the etna gem. These workflows are:

- Create Project: The initial command you run that lets you define multiple models and attributes.
- Add Model: Allows you to add a single model and its attributes to an existing project.
- Attribute Actions: A set of four different commands that let you add, rename, or update attributes on existing models.

Each workflow requires a JSON file input, and also offers a command to validate the structure of your JSON.

The general workflow for defining your project should look like:

1. Sketch out your project's model and attribute definitions in a spreadsheet or document.
2. Convert the models and attributes to the "create project" JSON format. You can follow the template if that would be helpful.
3. (optional) Iterate on the above and get feedback from the engineering team. Use the `validate` subcommands to check your JSON format.
4. Request that a project be created for you on Janus **and** Janus stage, with the same project names.
5. Run the "create project" command against our staging environment. You must be on the VPN for this step.
6. Visit Timur staging and verify the project structure appears correct (the Map tab may be particularly helpful for this). You must be on the VPN for this step.
7. Edit the project JSON as necessary.
8. Run the "create project" command against our production environment.
9. Visit Timur to verify the project structure appears correct.

### JSON Templates

Several JSON templates are available for your reference. They are in the GitHub repository at the URLs below, and you can also copy them out of the gem via a command like:

```bash
$ etna create_template create_project
A sample project JSON template has been provided at `project_template.json`.
```

When using the templates, make sure to remove all the comments, which are on lines starting with `//`. JSON format does not allow comments, and they are provided for initial explanation only.

| Action                | Template                                                                                                                                     | Etna Command                           |
| --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------- |
| Create Project Models | [create_project_template.json](https://github.com/mountetna/monoetna/blob/master/etna/lib/etna/templates/create_project_template.json)       | `$ etna create_template create_project`    |
| Add Single Model      | [add_model_template.json](https://github.com/mountetna/monoetna/blob/master/etna/lib/etna/templates/add_model_template.json)                 | `$ etna create_template add_model`         |
| Attribute Actions     | [attribute_actions_template.json](https://github.com/mountetna/monoetna/blob/master/etna/lib/etna/templates/attribute_actions_template.json) | `$ etna create_template attribute_actions` |
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
- desc
- display_name

The following values for `attribute_type` are supported in the JSON:

- boolean
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

The following validation `type` values are supported for each attribute that you want to define validation for:

- Array
- Range
- Regexp

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

#### Project

The following keys are required for the project definition:

- project_name
- project_name_full

### JSON Validation

The etna gem provides commands to validate the JSON structure for the three project definition workflows. They are:

```bash
$ etna administrate project validate __filepath__ --environment
$ etna administrate model validate __project_name__ __model_name__ __filepath__ --environment
$ etna administrate model attributes validate_actions __project_name__ __filepath__ --environment
```

These commands will output the results to the command line. Valid JSON files have a message like:

```bash
$ etna administrate project validate test_project.json
Project JSON is well-formatted!
```

Whereas invalid JSON files will report a list of errors, like:

```bash
$ etna administrate project validate missing_model_keys.json
Project JSON has 8 errors:
  * Parent model "" for document does not exist in project.
	Current models are ["assay_name", "document", "assay_pool", "patient", "project", "status", "symptom", "timepoint"].
  * Parent model "" for patient does not exist in project.
	Current models are ["assay_name", "document", "assay_pool", "patient", "project", "status", "symptom", "timepoint"].
  * Missing required key for model assay_name: "parent_link_type".
  * Missing required key for model document: "parent_model_name".
  * Missing required key for model document: "parent_link_type".
  * Missing required key for model assay_pool: "identifier".
  * Missing required key for model patient: "parent_model_name".
  * Missing required key for model project: "identifier".
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
