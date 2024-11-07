
## Building pipelines 

### Philosophy

- Create small generic jobs that can be re-used in different pipelines
- Try and leverage existing etna libs if possible
- Test these jobs with unit tests and not integration tests
- Let argo just act as the orchestrator/scheduler. 

### ETLJob super class

ETL jobs should subclass: `ETLJob` and implement the `pre`, `process`, and `post` methods.
This class has a whole bunch of useful methods .

### State management

It is often the case that your ETL job will need to manage some state. 

- If you need to store a file, in between pipeline runs, you should use a volume mount.

- For all other data you need to track, you should create a db table.

#### Argo 

Argo manages and collects metadata about the pipeline run, this includes things like start and end times, and status.
No need to create a table for this.

#### DB table

Your db table should always have an `argo_id` column. This is so we can match up the pipeline run with the db table.

Example columns for the cat ingestion pipeline:

- id
- argo_id
- last_scan
- num_files_to_update
- num_c4_files_updated
- num_metis_files_updated

Notice how these are all pipeline specific.

## Argo

### Submitting pipelines

### Viewing pipelines
