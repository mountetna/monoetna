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

### Config and run flow

- A user hits `save config` after they have made a param change OR they have written a ui file for a SINGLE UI STEP.
- A user cannot fill out multiple UI steps and then hit `save config`.
- They can keep hitting `save config` as many times as they want for a single UI STEP.

## Snakemake

### Let snakemake do its thing

Think of the Vulcan App as a wrapper around snakemake. We really want to let snakemake handle all the workflow logic, this includes things like caching, re-running jobs, etc.

### Config files

Snakemake workflows should use config files to specify parameters. This is the preferred way to run a snakemake workflow. Config files are globally scopped.

#### Default config file

- All workflows MUST have a default-config.json file. This file is called default-config.json. Snakemake will throw key errors if you do not provide a config file with all the requiered keys.
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

#### How target inference works

- We parse each snakemake file and extract the targets and their dependencies. So we would get an object that looks like:

```ruby
{"output/count_poem.txt"=>{"inputs"=>["output/poem.txt", "output/poem_2.txt"], "params"=>["count_bytes", "count_chars"]},
 "output/count_poem_2.txt"=>{"inputs"=>["output/poem.txt", "output/poem_2.txt"], "params"=>["count_bytes", "count_chars"]},
 "output/arithmetic.txt"=>{"inputs"=>["output/count_poem.txt", "output/count_poem_2.txt"], "params"=>["add", "add_and_multiply_by"]},
 "output/check.txt"=>{"inputs"=>["output/arithmetic.txt"], "params"=>[]},
 "output/summary.txt"=>{"inputs"=>["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt", "output/ui_job_one.txt", "output/ui_job_two.txt"], "params"=>[]},
 "output/final.txt"=>{"inputs"=>["output/ui_summary.txt"], "params"=>[]}}
 ```

- Now give a set of parameters and a set of available files, we can figure out which targets need to be run.

### Halting problem

The probem:

```Workflow:
params: [some gene expression files]

comp-step1: reading that file, and outputting a file of string[] (gene names)
 -input: params: [some gene expression files]
 -output: files: all-genes.txt

ui-step2: dropdown of those gene names
 -inputs: files: all-genes.txt
 -outputs: files: chosen-gene.txt

comp-step3: goes back to the file and outputs a file of two numbers: min(chosen-gene-expression), max(chosen-gene-expression)
 -inputs: files: chosen-gene.txt
 -outputs: files: gene-min-max.txt 

ui-step4: slider for user to pick a value between that min and max
 -inputs: files: gene-min-max.txt
 -outputs: files: chosen-cutoff.txt

comp-step5: Some other stuff
 -inputs: files: chosen-cutoff.txt
```

1) User had chosen geneX in ui-step2 (which produces and sends chosen-gene.txt, (then calls #save_config)), hits run (comp-step3 is now run, so produces gene-min-max.txt), then chosen some cutoff in ui-step4 (which produces and sends chosen-cutoff.txt, (then calls #save_config)).  Maybe they moved further, maybe not.
ouputs/
  chosen-gene.txt
  gene-min-max.txt
  chosen-cutoff.txt
2) But now has gone back to re-select geneY in ui-step2 (which produces and sends a DIFFERENT chosen-gene.txt, (then calls #save_config)).
outputs/
  chosen-gene.txt* <- snakemake notices this is different && knows comp-step3, ui-step4, comp-step5 would be set to be recreated
  gene-min-max.txt
  chosen-cutoff.txt <- BUT snakemake doesn't understand that this ui selection file output might be incompatible for geneY.  We've chosen that the path is to delete downstream ui outputs so that snakemake can nolonger blow right past this point.

(The potential for error if an old chosen-cutoff.txt is used: say geneX expression goes as high as 10000 tpm, but geneY expression only goes up to 100, and in 1) the cutoff chosen was 1500 that is out-of-range for geneY and likely leads to errors in comp-step5)

The solution:

- Delete the file outputs of scheduled/downstream ui-steps during or after this accounting & before a user has ability to trigger a next run

More details:

There are two cases we need to account for in order to properly remove UI files:

1. The params that have changed between the last config save and the current save, via a compute job.
2. The ui files that have just been written to the workspace, via a UI step.

The somewhat tricky part about this is we need to make sure we do not delete any upstream ui files that have just been written.

Ie: Lets say we've run this pipeline and generated files for each step:

ui_job_1 -> comp_job_2 -> ui_job_3

Now lets say we change a param in comp_job_2. We do not want to delete ui_job_1, since it is an upstream dependency of ui_job_2.
We just want to delete the files associated with ui_job_3.

In order to solve this problem, we build a file tree and use the two UI params `uiFilesSent` and `paramsChanged` to remove any downstream UI files.

# SSH

### Using the c4 env

- Specify username and password in the config.yaml

Ex:

```yaml
  :ssh: &default_ssh
    :host: ahost
    :username: user2
    :password: some-password
    :use_ssh_config: false
```

### Connecting to Prod

- Use a ssh config file and SSH key pair for the host you want to connect to

```yaml
  :ssh: &default_ssh
    :host: ahost
    :username: user2
    :password: 
    :use_ssh_config: true
```

### Conda

TODO:

```bashrc
# Define your conda paths at the top
CONDA_BASE_DIR="/software/c4/cbi/software/miniforge3-24.11.0-0"  # Base directory where miniforge is installed
CONDA_SH_PATH="$CONDA_BASE_DIR/etc/profile.d/conda.sh"

# Lazy load conda to improve shell startup time
conda() {
    # Remove the function after first execution
    unset -f conda
    
    # Load required modules first
    module load CBI miniforge3/24.11.0-0
    
    # Initialize conda
    . "$CONDA_SH_PATH"
    
    # Execute the original command
    conda "$@"
}
```



