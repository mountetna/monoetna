## Building pipelines 

### Data-eng Philosophy

- Create small generic jobs that can be re-used in different workflows 
- Each job should have a single responsibility
- Try and leverage existing etna libs
- Test these jobs with unit tests and not integration tests
- Let argo just act as the orchestrator/scheduler. 

### Directory structure 

- Build argo workflows in the `polyphemus/workflows` directory.
- Build etl jobs in the `polyphemus/lib/data_eng/` directory.

### Data model

- `configs` table: This table stores the workflow configs.
- `runs` table: This table stores the workflow runs.
- `run_metadata` table: This table stores the workflow run metadata.

A note on the configs table and ids.

There are two ids in the configs table.

1. `config_id`: This is the id of the config and does not change as multiple version of the configs are created.  Along with `version_number` it is the composite key for this table. 
2. `id`: This is the id of the row in the configs table. It is a unique identifier for a given workflow run.

`runs` and `run_metadata` tables both have a foreign key to the `configs` table via the `config_id` and `version_number` columns of the `configs` table.


### ETLJob super class

ETL jobs should subclass: `ETLJob` and implement the `pre`, `process`, and `post` methods.

### WorkflowManifest class

This class is responsible for definining:

- The name of the workflow
- The schema for the workflow config,
    - used for the UI to render components 
    - used for the definition and validation of the config
- The runtime params for the workflow
- The secrets for the workflow

It is intentionally separate from ETL jobs. Ideally the ETL jobs should be agnostic of the workflow they are a part of.


### RunJob command

We use the `RunJob` command to run an ETL job in the etna container.

bin/polyphemus run_job <workflow> <job_name> <other args>

### Language 

- Configs/Params: These are key value pairs that a *user* sets as the initial input for the ETL jobs. These can be configured in a UI by user in polyphemues, or by a simple json file that is then stored in db table. Things like metis directories, etc...

- Secrets: These are key value pairs that a *user* sets needed to run the ETL job. Things like ssh credentials, sftp credentials, etc.

- Workflow state: These are are values that the *workflow* and *jobs* generates for a given run. Things like how many files were updated, etc. Each workflow has it's own db table.

- Workflow meta-data: These are values that *argo* (or another orchestrator) generates for a given run. Things like how long the workflow took to run, failures, retries, etc. We can fetch this from argo. Argo generates a unique id for each run. 

### Polyphemus

#### ETL Jobs

ETL jobs use the polyphemus database for:

1. To retrieve configs for a workflow.
2. To manage pipeline state in the `runs` table.

The ETL jobs achieve this by using the Polyphemus client to make API calls.

#### Run metadata

There is a command called `WriteRunMetadata` that is used to write runtime metadata to the `run_metadata` table.
We use an argo onExit hook to run this once the workflow has completed.

### State management

#### Job state

It is often the case that your ETL job will need to manage some state. 
- If you need to store a file, in between pipeline runs, you should use a volume mount.
- For all other data you need to track, you should managed state in the `runs` table in the polyphemus database.

## Argo

### Submitting workflows

``` argo submit -f <path_to_workflow_yaml> -p config_id=<config_id> -p version_number=<version_number>```

### Viewing workflow runs
