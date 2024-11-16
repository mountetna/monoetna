## Building pipelines 

### Data-eng Philosophy

- Create small generic jobs that can be re-used in different pipelines
- Each job should have a single responsibility
- Each job should manage its own state in a database table
- Try and leverage existing etna
- Test these jobs with unit tests and not integration tests
- Let argo just act as the orchestrator/scheduler. 

### ETLJob super class

ETL jobs should subclass: `ETLJob` and implement the `pre`, `process`, and `post` methods.

### RunJob command

We use the `RunJob` command to run an ETL job in the etna container.

bin/etna run_job <pipeline> <job_name> <other args>

### Important concepts

- Configs: These are key value pairs that a *user* sets as the initial input for the ETL jobs. These can be configured in a UI by user in polyphemues, or by a simple json file that is then stored in the etl_configs table. Things like metis directories, etc...

- Secrets: These are key value pairs that a *user* sets needed to run the ETL job. Things like ssh credentials, sftp credentials, etc.

- Pipeline state: These are are values that the *pipeline* generates for a given pipeline run. Things like how many files were updated, etc.

- Pipeline meta-data: These are values that *argo* generates for a given pipeline run. Things like how long the pipeline took to run, failures, retries, etc.


### Relationship to Polyphemus

ETL jobs use the polyphemus database for 2 reasons:

1. To retrieve params for a pipeline from the `etl_configs` table.
2. To manage pipeline state in a specific db table.

The ETLJob itself is independant and unaware of polyphemus. Its only input is a secrets and a config. Which is the `config` and `secrets hash from the `etl_configs` table. And the ETL Job also doesn't know anything about the location of the 

### Data model


`etl_configs` <- 1 to many <- `pipeline_table`




#### ETL Configs Table

At the beginning of every job, we query the `polyphemus` db to retrieve the config for the pipeline. 

Some important things to note about the `etl_configs` table:

- `config` is a json blob of key/value pairs. This is how we allow each job to define it's own params.


The `etl_configs` table is where you manage the params for a pipeline. 


### Directory structure 

Since the etl jobs are just ruby classes, and leverage much of the etna lib, we can keep them in the `etna/lib` directory. 

- Build argo pipelines in the `etna/pipelines/argo` directory.
- Build etl jobs in the `etna/lib/jobs` directory.

### State management

#### Job state

It is often the case that your ETL job will need to manage some state. 
- If you need to store a file, in between pipeline runs, you should use a volume mount.
- For all other data you need to track, you should create a db table in the polyphemus database.

Ex: for the cat ingestion pipeline, we have a db table called `cat_ingestion_jobs`

Your db table should always have an `argo_id` column. This is so we can match up the pipeline run with the db table.

Example columns for the cat ingestion pipeline:

- id
- etl_configs_id
- argo_id
- last_scan
- num_files_to_update
- num_c4_files_updated
- num_metis_files_updated

Notice how these are all job specific. No need to collect meta-data about the pipeline run. Argo manages this.

## Argo

### Submitting pipelines

### Viewing pipelines
