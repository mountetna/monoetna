cwlVersion: v1.1
class: Workflow

inputs:
  1_Data_Source_magma_query__queryTerms:
    type: string
    label: "query input"
    doc: "The set of terms for constructing a magma query towards the data of interest.  Suggestion: Use the Query page of Timur to build the proper dataframe, and then simply copy over the 'query'-chunk presented in a yellow box there!"
  1_Data_Source_magma_query__user_columns:
    type: string
    label: "user_columns input"
    doc: "Optional array or comma-separated set of *quoted* strings for renaming columns of the returned data.  Suggestion: Use the Query page of Timur to build the proper dataframe, and then simply copy over the 'user_columns'-chunk presented there!"
  1_Data_Source_magma_query__expand_matrices:
    type: boolean
    label: "expand_matrices input"
    default: true
    doc: "Whether to expand matrix attributes into individual columns, with one matrix data point per column.  In most cases, you'll want to leave this checked."
  
outputs:
  the_plot:
    type: File
    outputSource: make_plot/plot.json

steps:
  get_data:
    run: scripts/VIZ_query_df.cwl
    label: 'Fetch Data'
    in:
      queryTerms: 1_Data_Source_magma_query__queryTerms
      user_columns: 1_Data_Source_magma_query__user_columns
      expand_matrices: 1_Data_Source_magma_query__expand_matrices
    out: [data_frame]
  fill_plot_options:
    run: ui-queries/y-plotly.cwl
    label: 'Set plot options'
    in:
      data_frame: get_data/data_frame
    out: [plot_setup]
  make_plot:
    run: scripts/make_yplot.cwl
    label: 'Create Y Plot'
    in:
      plot_setup: fill_plot_options/plot_setup
      data_frame: get_data/data_frame
    out: [plot.json]
  show_plot:
    run: ui-outputs/plotly.cwl
    in:
      a: make_plot/plot.json
    out: []
    label: 'Display Y Plot'
