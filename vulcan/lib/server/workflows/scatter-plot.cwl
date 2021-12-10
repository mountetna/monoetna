cwlVersion: v1.1
class: Workflow

inputs:
  1_Data_Source__queryTerms:
    type: string
    label: "queryTerms string"
    doc: "The set of terms for cunstructing a magma query towards the data of interest.  Suggestion: Use the Query page of Timur to build the proper dataframe, and then simply copy over the chunk presented there!"

outputs:
  the_plot:
    type: File
    outputSource: plot_scatter/plot.json

steps:
  get_data:
    run: scripts/VIZ_query_df.cwl
    label: 'Fetch Data'
    in:
      queryTerms: 1_Data_Source__queryTerms
    out: [data_frame]
  fill_plot_options:
    run: ui-queries/scatter-plotly.cwl
    label: 'Set plot options'
    in:
      data_frame: get_data/data_frame
    out: [plot_setup]
  plot_scatter:
    run: scripts/make_scatter.cwl
    label: 'Create Scatter Plot'
    in:
      plot_setup: fill_plot_options/plot_setup
      data_frame: get_data/data_frame
    out: [plot.json]
  show_scatter_plot:
    run: ui-outputs/plotly.cwl
    in:
      a: plot_scatter/plot.json
    out: []
    label: 'Display Scatter Plot'
