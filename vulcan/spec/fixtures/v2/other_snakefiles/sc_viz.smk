configfile: "config.yaml"

rule all:
    input:
        thumbnail="output/thumbnail.png",
        plot="output/plot.out"

rule get_dataset_and_summarize:
    params:
        the_dataset_name=config["dataset_name"]
    output:
        scdata="output/scdata.Rds",
        plotting_options="output/plotting_options.json",
        discrete_metadata_summary="output/discrete_metadata_summary.json",
        all_opts="output/all_opts.json",
        continuous_opts="output/continuous_opts.json",
        discrete_opts="output/discrete_opts.json",
        reduction_opts="output/reduction_opts.json"
    singularity:
        "/dscolab/vulcan/containers/archimedes-r.sif"
    script:
        "scripts/get_dataset_and_summarize.R"

# For UI elements that produce outputs used in a snakemake step, we just specify the input/output, so snakemake can infer the dag
rule plot_setup_ui:
    input:
        plotting_options="output/plotting_options.json"
    output:
        plot_settings="output/plot_setup.json"

rule make_plot:
    input:
        scdata="output/scdata.Rds",
        plot_setup="output/plot_setup.json",
        plotting_options="output/plotting_options.json"
    output:
        plot_out="output/plot.out",
        thumbnail="output/thumbnail.png",
        plot_Rds="output/plot.Rds"
    singularity:
        "/dscolab/vulcan/containers/archimedes-r.sif"
    script:
        "scripts/make_dittoSeq_plot.R"