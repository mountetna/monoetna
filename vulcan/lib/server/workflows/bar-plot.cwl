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
    outputSource: make_plot/plot.json

steps:
  get_data:
    run: scripts/VIZ_query_df.cwl
    label: 'Fetch Data'
    in:
      queryTerms: 1_Data_Source__queryTerms
    out: [data_frame]
  fill_plot_options:
    run: ui-queries/bar-plotly.cwl
    label: 'Set plot options'
    in:
      data_frame: get_data/data_frame
    out: [plot_setup]
  make_plot:
    run: scripts/make_barplot.cwl
    label: 'Create Bar Plot'
    in:
      plot_setup: fill_plot_options/plot_setup
      data_frame: get_data/data_frame
    out: [plot.json]
  show_plot:
    run: ui-outputs/plotly.cwl
    in:
      a: make_plot/plot.json
    out: []
    label: 'Display Bar Plot'
