cwlVersion: v1.1
class: Workflow

outputs:
  thumbnail:
    type: File
    format: image/png
    outputSource: make_umap/plot.png

steps:
  queryMagma:
    run: scripts/retrieve_sc_seq_datasets.cwl
    label: 'Fetch dataset options'
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
    out: [plot_options, sc_data, metadata_summary, continuous_metas, discrete_metas]
  plot_setup:
    run: ui-queries/single-cell-plot-any.cwl
    label: 'Set plot options'
    doc: "Options here determine both what type of plot to make, and how to set that plot up. For addtional details, see the 'Visualization with Vulcan' section of the Vulcan's documentation, acccessible via the 'Help' button at the top of this page. This particular instance of the Plot Configuration Interface constitutes a version with preset values for plot-type (scatter_plot), X-Axis Data (UMAP_1), Y-Axis Data (UMAP_2), and Color Data (chosen above)."
    in:
      plot_options: get_dataset_and_summarize/plot_options
      metadata_summary: get_dataset_and_summarize/metadata_summary
      continuous_metas: get_dataset_and_summarize/continuous_metas
      discrete_metas: get_dataset_and_summarize/discrete_metas
    out: [plot_setup]
  make_umap:
    run: scripts/make_plot.cwl
    label: 'Create UMAP plot'
    in:
      plot_setup: user_plot_setup/plot_setup
      data_frame: prep_umap_plot_data/data_frame
    out: [plot.json, plot.png]
  show_umap_plot:
    run: ui-outputs/plotly.cwl
    label: 'Display UMAP'
    in:
      a: make_umap/plot.json
    out: []
