https://github.com/mountetna/vulcan/workflows/Run%20JavaScript%20tests/badge.svg https://github.com/mountetna/vulcan/workflows/Run%20Ruby%20tests/badge.svg

# README

Vulcan is a workflow engine for running bioinformatic data pipelines written in Snakemake.

## Vulcan config files

### Snakemake jobs

### UI jobs

# Backend 

## Backend TLDR

- We use the cancer center HPC and we leverage it for running workflows.
- We have a networked disk on the HPC that is used for creating workspaces.
- A workspace is just a directory (ie: /app/vulcan/workspaces/123/my-workspace) where a snakemake workflow is copied and run.
- We invoke ssh commands from the Vulcan app to the HPC to perform all actions.

## Snakemake

### Let snakemake do its thing

Think of the Vulcan App as a wrapper around snakemake. We really want to let snakemake handle all the workflow logic, this includes things like caching, re-running jobs, etc.

### Config files

Snakemake workflows should use config files to specify parameters. This is the preferred way to run a snakemake workflow. Config files are globally scopped.

#### Default config file

- All workflows MUST have a default-config.yaml This file is called default-config.yaml. Snakemake will throw key errors if you do not provide a config file with all the requiered keys.
- This file is also used to generate a dag.
- Snakemake workflows should not include the directive, `configfile`, whenever you run a snakemake workflow, always use the `--config-file` flag. This makes it more explicity what config file is being used.

#### Config values in the param directive

- It is important that all config values must be passed in the param directive. This is both more explicit and allows snakemake to figure out which targets need to be run and re-run. See `--rerun-triggers` for more details.   

###  Dependency resolution

Recall that in bioinformatics, files are the first class citizens. As such snakemake was written and operates around files being inputs and outputs of jobs. Files that need to be genereted are called targets

There is a lot of code written in the backend around "Target Inference". The goal of this code is to determine which targets (files) in a snakemake workflow should be run based on the provided parameters and available files.

In general snakemake handles dependency resolution by itself, that is, given inputs/outputs/params, snakemake will determine which targets need to be run.

We do allow snakemake to figure out which targets need to be run or re-run, however we restrict the target files to a specific subset of files. We do this for the following reasons:

- If a user only fills out params for two rules, and then wants to inspect the output of those rules, like via a Vulcan UI job, then we only want to run those two rules (and their corresponding targets) and stop. If we do not restrict to these targets, and since we the entire config is defined before hand (hence snakemake thinks it has all the necessary config values), snakemake may run more than just the two rules.

- Ensures deterministic behavior.

- Faster, we only run the targets we need to run.

NOTE: It is easier and more efficient for snakemake if we tell snakemake which targets to run instead of using rules.

