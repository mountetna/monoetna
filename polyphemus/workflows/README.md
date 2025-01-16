## Building pipelines 

### Data-eng Philosophy

- Create small generic jobs that can be re-used in different workflows 
- Each job should have a single responsibility
- Try and leverage existing etna libs
- Test these jobs with unit tests and not integration tests
- Let argo just act as the orchestrator/scheduler. 

### Directory structure 

- Build argo workflows in the `polyphemus/workflows` directory.
- Build etl jobs, related clients, and workflow manifests in the `polyphemus/lib/data_eng/` directory.

### Language 

- Configs/Params: These are key value pairs that a *user* sets as the initial input for the ETL jobs. These can be configured in a UI by user in polyphemues, or by a simple json file that is then stored in db table. Things like metis directories, etc...

- Secrets: These are key value pairs that a *user* sets needed to run the ETL job. Things like ssh credentials, sftp credentials, etc.

- Workflow state: These are are values that the *workflow* and *jobs* generates for a given run. Things like how many files were updated, etc. Each workflow has it's own db table.

- Workflow run meta-data: These are values that *argo* (or another orchestrator) generates for a given run. Things like how long the workflow took to run, failures, retries, etc. We can fetch this from argo. Argo generates a unique id for each run. 

- Workflow manifests: These are Ruby classes that encapsulate some meta-data about the workflow.

This class is responsible for definining:

- The name of the workflow
- The schema for the workflow config,
    - used for the UI to render components 
    - used for the definition and validation of the config
- The runtime params for the workflow
- The secrets for the workflow

Note they are intentionally separate from ETL jobs. The ETL jobs should be agnostic of the workflow they are a part of.

### Data model

- `configs` table: This table stores the workflow configs.
- `runs` table: This table stores the workflow runs. 
- `runtime_config` table: This table stores the runtime config for a workflow.

A note on the configs table and ids.

There are two ids in the configs table.

1. `config_id`: This is the id of the config and does not change as multiple version of the configs are created.  Along with `version_number` it is the composite key for this table. 
2. `id`: This is the id of the row in the configs table. It is a unique identifier for a given workflow run.

`runs` has a foreign key to the `configs` table via the `config_id` and `version_number` columns of the `configs` table. We could have also used the `id` column of the `configs` table, but this is less explicit, and it is forces us to think of a "run" as a single instantiation of config and version number.

It is also worth noting that the `runtime_config` table should only ever have one row per config. That row is always updated.

#### Workflow names and UID

Argo generates two useful idenfiers when a workflow is run:

- workflow.name: This is a human readable name and is unique for the duration of the worfklow. It MAY be re-used after the workflow has been cleared.

- workflow.uid: Th unique Kubernetes-assigned UID

The workflow.uid is known as the `run_id` in our db. 

### ETLJob super class

ETL jobs should subclass: `ETLJob` and implement the `pre`, `process`, and `post` methods.
You care store state in the context hash between methods.

### RunJob command

We use the `RunJob` command to invoke ETLJobs.

bin/polyphemus run_job <workflow> <job_name> <config_id> <version_number>

Note: If you want to invoke the `SFTPFileDiscoveryJob`, in the workflow yaml, you must use the string `sftp_file_discovery`.
Ex: 


### ETL Jobs

ETL jobs use the polyphemus database for:

1. To retrieve configs for a workflow from the `configs` table.
2. To manage pipeline state in the `runs` table.

The ETL jobs achieve this by using the Polyphemus client to make API calls.

### Run metadata

There is a command called `RecordOrchestratorMetadata` that is used to write runtime metadata to the `runs` table.
We use an argo onExit hook to run this once the workflow has completed or exited for any reason.

### Run intervals

There is a column in the `runtime_configs` table called `run_interval`. This is used to determine how often the workflow should be run.

WIP -> We can probably write some code and use CronWorkflows to manage this.

### State management

#### Job state

It is often the case that your ETL job will need to manage some state. 
- If you need to store a file, in between pipeline runs, you should use a volume mount.
- For all other data you need to track, you should managed state in the `runs` table under the `state` column in the polyphemus database.

## Argo

### Submitting workflows

``` argo submit -f <path_to_workflow_yaml> -p config_id=<config_id> -p version_number=<version_number>```

### Viewing workflow runs
