cwlVersion: v1.1
class: Workflow

inputs:
  1_Object_Selection__dataset_name:
    type: string
    default: ''
    label: 'name of record to target'
    doc: 'Provide the identifier ("name" attribute) of the sc_seq_dataset record you wish to explore here. Note that we definitely plan to improve this particular selection method in the future!'
outputs:
  thumbnail:
    type: File
    format: image/png
    outputSource: make_plot/plot.png

steps:
  get_dataset_and_summarize:
    run: scripts/get_dataset_and_summarize.cwl
    label: 'Retrieve dataset & Determine plotting options'
    in:
      dataset_name: 1_Object_Selection__dataset_name
    out: [plotting_options, scdata, discrete_metadata_summary, all_opts, continuous_opts, discrete_opts, reduction_opts]
  plot_setup:
    run: ui-queries/any-dittoseq.cwl
    label: 'Set plot options'
    doc: "Options here determine both what type of plot to make, and how to set that plot up. For addtional details on individual inputs, see the 'Inputs of the Plot Configuration Interface' section of Vulcan's 'Help'-page documentation OR dittoSeq's own documentation. Any inputs without an exact label match within that Vulcan 'Help'-page table will map directly to some dittoSeq input for the function with the same name as your chosen plot type. You can find dittoSeq's documentation from either within R itself, run `?dittoSeq::<visualization-name>`, or from the RDRR.io page that pops up when you google the package."
    in:
      data_frame: get_dataset_and_summarize/discrete_metadata_summary
      continuous_cols: get_dataset_and_summarize/continuous_opts
      discrete_cols: get_dataset_and_summarize/discrete_opts
      all_cols: get_dataset_and_summarize/all_opts
      reduction_opts: get_dataset_and_summarize/reduction_opts
    out: [plot_setup]
  make_plot:
    run: scripts/make_dittoSeq_plot.cwl
    label: 'Create Plot'
    in:
      plot_setup: plot_setup/plot_setup
      scdata: get_dataset_and_summarize/scdata
      plotting_options: get_dataset_and_summarize/plotting_options
    out: [plot.out, legend.png, plot.png, plot.Rds]
  download_plot:
    run: ui-outputs/link.cwl
    in:
      a: make_plot/plot.Rds
    out: []
    label: 'Download plot object'
  show_plot:
    run: ui-outputs/plot.cwl
    label: 'Display Plot'
    in:
      a: make_plot/plot.out
    out: []


