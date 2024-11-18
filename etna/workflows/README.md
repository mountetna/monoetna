## Building pipelines 

### Data-eng Philosophy

- Create small generic jobs that can be re-used in different workflows 
- Each job should have a single responsibility
- Try and leverage existing etna libs
- Test these jobs with unit tests and not integration tests
- Let argo just act as the orchestrator/scheduler. 

### Directory structure 

Since the etl jobs are just ruby classes, and leverage much of the etna lib, we can keep them in the `etna/lib` directory. 

- Build argo pipelines in the `etna/pipelines/argo` directory.
- Build etl jobs in the `etna/lib/jobs` directory.

### ETLJob super class

ETL jobs should subclass: `ETLJob` and implement the `pre`, `process`, and `post` methods.

### RunJob command

We use the `RunJob` command to run an ETL job in the etna container.

bin/etna run_job <workflow> <job_name> <other args>

### Language 

- Configs/Params: These are key value pairs that a *user* sets as the initial input for the ETL jobs. These can be configured in a UI by user in polyphemues, or by a simple json file that is then stored in db table. Things like metis directories, etc...

- Secrets: These are key value pairs that a *user* sets needed to run the ETL job. Things like ssh credentials, sftp credentials, etc.

- Workflow state: These are are values that the *workflow* and *jobs* generates for a given run. Things like how many files were updated, etc. Each workflow has it's own db table.

- Workflow meta-data: These are values that *argo* (or another orchestrator) generates for a given run. Things like how long the workflow took to run, failures, retries, etc. We can fetch this from argo. Argo generates a unique id for each run. 

### Relationship to Polyphemus

ETL jobs should use the polyphemus database for 2 reasons:

1. To retrieve configs/params for a pipeline from a specific table.
2. To manage pipeline state in a specific db table.

The ETL jobs achieve this by using the Polyphemus client and the following api calls:

```      
def get_workflow(project_name, workflow_name, revision: "latest")
end

def get_workflow_state(argo_id)
end

def update_workflow_state(argo_id, state)
end
```

### Existing data model

Currently configs/params for a pipeline are stored in the `etl_configs.configs` table. 

### State management

#### Job state

It is often the case that your ETL job will need to manage some state. 
- If you need to store a file, in between pipeline runs, you should use a volume mount.
- For all other data you need to track, you should create a db table in the polyphemus database.

Ex: for the cat ingestion workflow, we have a db table called `workflow_cat_ingestion`

Your db table should always have:
- an `argo_id` column. 
- a `workflow_id` column. 

Example columns for the cat ingestion pipeline:

- id
- workflow_id
- argo_id
- last_scan
- num_files_to_update
- num_c4_files_updated
- num_metis_files_updated

Notice how these are all job specific. No need to collect meta-data about the pipeline run. Argo manages this.

## Argo

### Submitting pipelines

### Viewing pipelines
