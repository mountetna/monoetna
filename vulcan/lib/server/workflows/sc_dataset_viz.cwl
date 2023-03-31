cwlVersion: v1.1
class: Workflow

inputs:
  1_Cell_Filtering__min_nCounts:
    type: int
    default: 200
    label: 'minimum UMIs'
    doc: 'Minimum number of unique transcript reads per cell. Cells with fewer will be discarded. Recommended range: 200 (lower quality data, macrophage retention) to 1000 (very high quality data). Set to 0 to skip.'

outputs:
  thumbnail:
    type: File
    format: image/png
    outputSource: make_plot/plot.png

steps:
  queryMagma:
    run: scripts/retrieve_sc_seq_datasets.cwl
    label: 'Fetch dataset options'
    in:
      first: 1_Cell_Filtering__min_nCounts
    out: [dataset_options]
  selectDataset:
    run: ui-queries/select-autocomplete.cwl
    label: 'Dataset selection'
    in:
      a: queryMagma/dataset_options
    out: [selected_record]
  get_dataset_and_summarize:
    run: scripts/get_dataset_and_summarize.cwl
    label: 'Retrieve dataset & Determine plotting options'
    in:
      dataset_name: selectDataset/selected_record
    out: [plotting_options, scdata, discrete_metadata_summary, all_opts, continuous_opts, discrete_opts, reduction_opts]
  plot_setup:
    run: ui-queries/any-dittoseq.cwl
    label: 'Set plot options'
    doc: "Options here determine both what type of plot to make, and how to set that plot up. For addtional details, see the 'Visualization with Vulcan' section of the Vulcan's documentation, acccessible via the 'Help' button at the top of this page. This particular instance of the Plot Configuration Interface constitutes a version with preset values for plot-type (scatter_plot), X-Axis Data (UMAP_1), Y-Axis Data (UMAP_2), and Color Data (chosen above)."
    in:
      data_frame: get_dataset_and_summarize/discrete_metadata_summary
      continuous_cols: get_dataset_and_summarize/continuous_opts
      discrete_cols: get_dataset_and_summarize/discrete_opts
      all_cols: get_dataset_and_summarize/all_opts
      reduction_opts: get_dataset_and_summarize/reduction_opts
    out: [plot_setup]
  make_plot:
    run: scripts/make_plot.cwl
    label: 'Create Plot'
    in:
      plot_setup: plot_setup/plot_setup
      scdata: get_dataset_and_summarize/scdata
      plotting_options: get_dataset_and_summarize/plotting_options
    out: [plot.json, plot.png]
  show_plot:
    run: ui-outputs/plotly.cwl
    label: 'Display Plot'
    in:
      a: make_plot/plot.json
    out: []


