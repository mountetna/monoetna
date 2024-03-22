import subprocess
import os
import argparse
import yaml
import json


def run_snakemake(local, profile, workflow_profile, singularity_args):

    snakemake_cmd = f"snakemake"

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

def parse_config():
    config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.yaml')
    with open(config_path, 'r') as file:
        config_data = yaml.safe_load(file)
    print(json.dumps(config_data, indent=2))

def main():
    parser = argparse.ArgumentParser(description='Run Snakemake pipelines')
    parser.add_argument('--run', action='store_true', help='Run pipelines')
    parser.add_argument('--local', action='store_true', help='Do not use Slurm for job scheduling')
    parser.add_argument('--profile', default='profiles/generic/', type=str, help='Full path to the Snakemake profile directory')
    parser.add_argument('--workflow-profile', default='', type=str, help='Full path to workflow specific profiles - must be directory.')
    parser.add_argument('--singularity-args', default='', type=str, help='Additional arguments to pass to Singularity')
    parser.add_argument('--parse-config', action='store_true', help='Parse a fixed YAML config file and print it as JSON')

    args = parser.parse_args()

    if args.parse_config:
        parse_config()
    elif args.run:
        run_snakemake(args.pipeline, args.local, args.profile, args.workflow_profile, args.singularity_args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
