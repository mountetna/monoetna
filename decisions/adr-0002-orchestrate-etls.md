# ADR-0002: Orchestrate ETLs

## Status

Accepted

## Context

We want an ETL orchestrator to manage and schedule tasks such as:

- ingestion from a data source (i.e. Box, CAT) to the data library.
- automate linking of data from Metis files to Timur.
- backup databases.

Ideally, we imagine this orchestrator may be also able to help with seemingly similar tasks:

- run analysis pipelines on C4
- run on-demand computation from Vulcan (possibly on C4, possibly not)

We would also like to have an interface for data scientists and other users to track progress of ETL tasks, define their own ETLs, and start ETL runs.

Currently, we use a bespoke ETL and cursor architecture written in Ruby and run via Polyphemus. Some of the challenges around this are:

- Only engineers can write new ETLs, since they are written with Ruby.
- It is not schedulable outside of Polyphemus (REDCap ingestion) or `systemctl` jobs. Scheduling changes require editing Chef configurations.
- Deploying new ETLs requires a code deployment. This can lead to long feedback loops.
- Difficult to test and develop outside of a production environment, often requiring the author to manually create local test data.

## Decision

We decided to use Apache Airflow, since it is widely supported as an ETL orchestrator and has many plugins. It also includes a dashboard where users can track ETL progress, view logs, and trigger ETL runs. We will use `v2`, which has a REST API that we could use to integrate with other systems, like Vulcan.

We run two types of tasks in Airflow:

- Native, Python DAGs and tasks. These use a LocalExecutor, for simplicity, and could be swapped out later. These are generally newly written in Python, or ported over from Ruby ETLs.
- Containerized ETLs running Ruby code. These come in several flavors:

  - Polyphemus ETLs, that we have not ported to Python.
  - App-related tasks, like backing up an application database.

For the containerized ETLs, we use a custom Docker Swarm operator for these, not the included Airflow Docker operator. This is so that we can copy our application service configs for the app containers (i.e. networks, resource limitations, constraints, etc.), to make DAG task definitions simpler. We allow overwrites for some of the service configurations, since DAG tasks are not exact clones of the application definitions. For example, the `sync_GNE` task uses a Metis app container, but it requires more CPU and memory than the application itself, so in the Airflow DAG definition, we have increased the memory and CPU limits.

Also, the containerized ETLs are started with a service name like `DockerServiceOperator_<hash>`. In Portainer or the Docker Swarm CLI, you can then find these service tasks when they are running and check things like resource usage or logs.

### Restarting ETLs

Each ETL has a version number associated with it. To re-start an entire ETL, you can bump the version number in the DAG decorator, using the Airflow code editor. This causes a new DAG to be created.

- Note that for the containerized Polyphemus ETLs that use our Ruby ETL + Cursor framework (i.e. sync_GNE), restarting those requires going into a Polyphemus console (via Portainer, a `polyphemus_app` service) and running `$ ./bin/polyphemus etl <tab_to_get_name_of_etl> reset`.

### Project-scope

ETLs are project-scoped, which means users should only be able to see the DAGs associated with their Janus token in the main Airflow dashboard. This is controlled via a custom Authentication backend (`etna_cookie_auth.py`). However, we have not implemented similar controls in `airflow_code_editor`, so any user can edit Python DAGs for any project.

Project DAGs are organized in project folders for convenience.

You can create a new project folder through the Admin -> New Project button, or by saving a DAG file and putting the new folder name into the "save" path.

### DAG Backups

The DAGs are being backed up nightly to Metis via a `backup` task, so we have not enabled the remote GitHub push for the code editor.

### Vulcan + C4

For now, we have left out integration with Vulcan and C4. The requirements for this integration were not clear, and it's not clear that Vulcan integrating with Airflow would satisfy all the system requirements. For example, many Vulcan tasks are on-offs, whereas Airflow is designed to repeatedly run the same tasks over different time windows. We hope with more concrete use cases, we can determine what the right integration with Vulcan should be or if we need another orchestrator / scheduler for that type of work.

### Airflow Pool Configuration

Airflow does not do a resource check for available resources before launching tasks, including checking the status of Docker Swarm (though it also appears that Docker Swarm itself does not check node resource availability before starting services). So we've set the Airflow pool size to 10 to limit the number of simultaneously running tasks. This seems to work well for our Swarm size, letting long-running tasks complete without blocking shorter-lived tasks (have also tried a pool size of 5, which led to long-running tasks blocking all other tasks).

## Consequences

- Requiring that Airflow has the ability to launch Docker Swarm services means that the Airflow scheduler has to run on a Swarm manager node, so that it can create and run Swarm services. Also, since we are using a LocalExecutor for non-Docker tasks, they will run on the same node as the Airflow schedule. To accomodate this local workload plus the need for the scheduler to be on a manager node, we have designated `dsco1` as a beefy manager node dedicated solely to Airflow. I think if we used a Celery executor instead of the LocalExecutor, we could then isolate the Airflow scheduler service onto a lighter-weight manager node, and have Airflow Celery consumers exist on the beefier, worker nodes. This might also require our manager nodes to have more resources though, since Airflow itself is resource intensive. We have observed that even slightly overloading the manager nodes often leads to frozen / hung docker daemons, and it's possible that Airflow might overwhelm our current manager nodes' resources.

- Since containerized ETLs run as Docker Swarm services (i.e. `DockerServiceOperator_<hash>`), they are ephemeral, and it can be very difficult to get log information for short-lived / crashing / buggy tasks through Portainer, since they appear in the list of services for such a short period of time. It would be nice if those logs could persist somewhere, or if the log files were more easily findable even after a service stops running.

- ETLs should run in time-batches, per the basic Airflow design (with start_time and end_time windows). This means that some of our ETLs, which do not have good time information (i.e. Box and CAT ingestion ETLs), keep their own cursor information in an Airflow variable as an alternative to time batching. For example, the Box ingest ETL uses file size plus name to create a unique hash, which is stored in a variable and compared against a file list to determine which files have already been ingested. At some point we may run into size limitations of these variable blobs. There also seem to be concurrency issues where parallel tasks (i.e. in CAT ingestion) try to update the same variable at the same time, and one will get locked out. Task retries automatically overcome those collisions, but it would be nice to find a more robust solution.

- Because ETLs are project-scoped, we need to configure a new Airflow connection for each target project, with a specific task token created for each connection. This is additional configuration overhead and may lead to missed ETL runs when the tokens expire.

- Many ETLs take advantage of shared Python code, which then requires ETL authors to be familiar with the structure and methods available in the helper classes. For example, there is a set of methods for single-cell file linking. However, DAG authors really need to be familiar with the methods, their arguments, and the defaults, in order to create single-cell tasks in Airflow. We do provide documentation for the `Etna` provider, but it still seems like a barrier for the data scientists to authoring their own ETLs. Perhaps some templates to copy from (i.e. in `coprojects_template`) would help adoption.

- We currently have limited notifications or error reporting from Airflow to Slack or e-mail, only some success messages for ingest ETLs. We decided to minimize this type of noise, so that folks would not tune out Airflow notifications. However, when ETLs are stuck (because a Swarm node hangs), we are never made aware. Better monitoring and notifications could help here. Right now this requires manual checking of the Airflow UI or someone complaining about missing ingest data, before we know that ETLs are failing.

- Better logging and monitoring of ETLs would be helpful. It is hard to search through the existing task logs, especially if they run multiple times before they can complete. It's also unclear how to best monitor and log tasks that might run on C4 or from Vulcan.

- Many ETL tasks are difficult to test without production data, so we still have the issue of long feedback loops and difficulty testing. We now tend to push Airflow images directly to production to test, especially if production data is required. There may be better ways to handle these challenges.
