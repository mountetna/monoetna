# ADR-0007: Vulcan Workflows

## Status

Proposed

## Context

Two broad categories of work exist surrounding the data library.

The first are data processing pipelines (DPP), which begin with raw data and produce research-ready data, such as a set of mutation calls or a cell x gene counts matrix. These are usually extremely compute-intensive and are expected to run more-or-less once.

The second are data analysis pipelines (DAP), which operate on the research-ready data produced by the DPPs. They may or may not be compute-intensive, but are often run iteratively to refine analytic results or explore data. DAPs may be expected to produce images, tables and other artifacts as ways to introspect on data.

Currently, DPPs are created and maintained outside of the library, and DAPs are created and maintained as part of the Vulcan code base. In addition, DPPs are configured and run manually by data scientists on external infrastructure, while DAPs run within Vulcan.

While these two kinds of pipeline thus appear to be separate domains at the moment, there are also a variety of factors drawing them together. In both cases, the pipeline consists of a directed acyclic graph (DAG) of tasks; steps in the pipeline operate on files and yield files as output; each step ought to run in a container and must have the appropriate libs, packages, etc., available to it; the DAG must be run by a resource-aware scheduler that is able to queue work when constrained.

In other words, in terms of execution, the two kinds of pipeline have virtually identical requirements; it cannot be assumed that they should therefore be executed in different computing environments. Both DPPs and DAPs are capable of generating system requirements that could benefit from access to HPC computing resources.

## Decision

We wish to unify the space of workflows in the data library, and integrate the normal work of data scientists in our group (and others) into the library. To do so, we will introduce a new Workflow data store in Vulcan. The data store will consist of two parts, a database model and a git tree. Workflow records will associate at the project level (i.e., each project will have a list of workflows). There will be no workflow author; the workflow belongs to the project.

### Project workflows
The Vulcan API will allow creating and listing project workflows. Each workflow can be also be browsed via the Vulcan UI. Workflows can still be run to produce Figures in Vulcan, which allow the workflow to be parameterized and run via the Vulcan UI. The Vulcan "steps" UI will be improved to monitor progress of a running workflow and show logging.

The Vulcan API will also expose a git remote endpoint for each project workflow. When the workflow is created, it will be initialized to an empty tree or forked from an existing workflow. This is a normal git remote and can be pushed/pulled from as desired; the endpoint can possibly validate the push using an update hook.

### Using Snakemake for workflows

The format of a Vulcan workflow is a Snakemake workflow. This may include steps in python, R, etc. Workflows may define containers via a Snakefile and conda yaml. Vulcan will keep track of containers and re-use them across workflows.

Inputs to the workflow can come from two places: 1) the data library, via https get, and 2) the cache/ folder associated with a Vulcan figure. Similarly, outputs can be: 1) writes to the data library (Magma or Metis) and 2) writes to the cache/ directory of a given folder

Workflows will be executed by running Snakemake in the appropriate environment. Each job will contain a Snakefile, a config.yml, Vulcan inputs and an Etna token, and will generate output which will be copied into Vulcan's cache when the workflow completes or halts. Here we remain agnostic about how jobs will be dispatched.

### UI Queries

Currently the Vulcan orchestrator has special hooks for "ui query" steps. Since we cannot expect the Snakemake orchestator to understand this logic, we will instead separate these out into a new set of Vulcan endpoints, which allow a workflow to set a requirement for a UI query associated with a particular step in the workflow. To make use of this interface, a step in a Snakemake workflow which expects to output a UI query JSON (i.e., the JSON response from the user) will use a vulcan client interface to set a query requirement and exit, halting the workflow due to the unsatisfied output requirement. The Vulcan UI will then be able to determine an unresolved query and will show it to the user in the UI, making use of the existing query interfaces to resolve the query. This results in the response being written to the appropriate output location, and the Snakemake workflow resumes with the requirement satisfied.

### Using Vulcan via the Command-Line

Although it may eventually be attractive to run a given workflow in the Vulcan UI, the development loop in the UI is difficult. Files cannot be inspected directly, and other utilities, REPLs, etc, may not be accessible. Although the use of Snakemake workflows makes it easier to transfer Vulcan workflows to other environments and run them there, many dependencies (containers and inputs) may be more difficult to achieve than a git clone. To allow easier development, Vulcan will include a command-line client, which will run the workflow's Snakefile from the command-line in a given environment, allowing users to use Vulcan containers and attach inputs in local folders. Similarly the command-line client will allow a user to resolve a UI query step in the shell by posting the appropriate output JSON to the correct location.

### Launching REPLs

A Vulcan figure collecting a set of data outputs (which may include data files, tables, PNGs or interactive graphs) is a "final artifact" which is shared through the library. However, often research involves taking data and playing with it through a process of trial and error. In these cases it may be useful to be able to attach the figure outputs into a REPL environment. To facilitate this, each Vulcan figure will provide a button which will launch a Jupyter notebook or RStudio environment with the figure data attached, within a container appropriate to playing with the data (defined within the workflow).

## Consequences

The results of our work should:

1. Promote sharing of workflows via the library
2. Allow easy movement between a command-line environment and Vulcan
3. Allow the use of the Vulcan user interface for parameterizing interactive workflows.
4. Allow the easy sharing of data products of workflows via the library
5. Play easily with workflow output in a REPL

As a result of this work, we hope that data scientists will be encouraged to do experimental work using the library. This will produce a steady accumulation of intermediate results (figures containing plots and other data) as a project unfolds. It will also aid in the progressive refinement of research by making work methodology available alongside results, so others may not only learn method but even copy and run working code.
