import subprocess
import os
import argparse
import json

def transform_config(config):
    flat_config = []

    def flatten(d, parent_key=''):
        for k, v in d.items():
            new_key = f"{parent_key}.{k}" if parent_key else k
            if isinstance(v, dict):
                flatten(v, new_key)
            else:
                flat_config.append(f"{new_key}={v}")

    flatten(config)
    return ' '.join(flat_config)

def run_snakemake(config, until, local, profile, workflow_profile, singularity_args):
    snakemake_cmd = f"snakemake"

    if until:
        snakemake_cmd += f" --until {until}"

    if config:
        json_dict = json.loads(config)
        flat_config = transform_config(json_dict)
        snakemake_cmd += f" --config {flat_config}"

    # Use Slurm unless --local is specified
    if local:
        snakemake_cmd += " --cores 1"
    else:
        snakemake_cmd += f" --profile {os.path.abspath(profile)}"
        if workflow_profile:
            snakemake_cmd += f" --workflow-profile {os.path.abspath(workflow_profile)}"

    # Pass in custom args
    if singularity_args:
        snakemake_cmd += f" --use-singularity --singularity-args '{singularity_args}'"

    print(f"Running command: {snakemake_cmd}")
    subprocess.run(snakemake_cmd, shell=True)

def main():
    parser = argparse.ArgumentParser(description='Run Snakemake pipelines')
    parser.add_argument('--run', action='store_true', help='Run pipelines')
    parser.add_argument('--config', default='', type=str, help='Pass config arguments to snakemake. Must be json')
    parser.add_argument('--until', default='', type=str, help='Run snakemake until a specific job')
    parser.add_argument('--local', action='store_true', help='Do not use Slurm for job scheduling')
    parser.add_argument('--profile', default='profiles/generic/', type=str, help='Full path to the Snakemake profile directory')
    parser.add_argument('--workflow-profile', default='', type=str, help='Full path to workflow specific profiles - must be directory.')
    parser.add_argument('--singularity-args', default='', type=str, help='Additional arguments to pass to Singularity')

    args = parser.parse_args()

    if args.run:
        run_snakemake(args.config, args.until, args.local, args.profile, args.workflow_profile, args.singularity_args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
