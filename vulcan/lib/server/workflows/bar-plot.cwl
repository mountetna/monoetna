cwlVersion: v1.1
class: Workflow

inputs:
  1_Data_Source_magma_query__queryTerms:
    type: string
    label: "query input"
    doc: "The set of terms for constructing a magma query towards the data of interest.  Suggestion: Use the Query page of Timur to build the proper dataframe, and then simply copy over the 'query'-chunk presented in a green box there!"
  1_Data_Source_magma_query__user_columns:
    type: string
    label: "user_columns input"
    doc: "Optional comma-separated set of *quoted* strings for renaming columns of the returned data.  Suggestion: Use the Query page of Timur to build the proper dataframe, and then simply copy over the 'user_columns'-chunk presented there!"
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
    run: ui-queries/bar-plotly.cwl
    label: 'Set plot options'
    doc: "Selections here adjust the plot that will be created. 'X-Axis Data' and 'Y-Axis Data' dropdowns present column names of your chosen data as their options. Note on 'make': Any text left as, or selection of, 'make' refers the visualization machinery to fill in the associated feature with a default value. For example, for 'X-Axis Title', leaving this text input as 'make' will let this title be filled in with its default, which we've defined to be whatever you selected from the 'X-Axis Data' dropdown."
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
