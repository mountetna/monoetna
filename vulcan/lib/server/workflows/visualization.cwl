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
  thumbnail:
    type: File
    format: image/png
    outputSource: make_plot/plot.png

steps:
  get_data:
    run: scripts/VIZ_query_df.cwl
    label: 'Fetch Data'
    in:
      useQuery: 1_Data_Source_magma_query__useDataFromQuery
      queryTerms: 1_Data_Source_magma_query__queryTerms
      user_columns: 1_Data_Source_magma_query__user_columns
      expand_matrices: 1_Data_Source_magma_query__expand_matrices
    out: [data_frame]
  review_data:
    run: ui-queries/data-transformation.cwl
    label: 'Review your data frame'
    doc: "Review and edit your data frame as needed. Right click for an interactions menu where you can add/remove columns. Start a cell with '=' to create functions in an Excel-like manner."
    in:
      data_frame: get_data/data_frame
    out: [formulaic_data, calculated_data]
  extend_user_formulas:
    run: scripts/calc_data_frame.cwl
    label: 'Extend user formulas'
    in:
      original_data.json: get_data/data_frame
      user_data.json: review_data/formulaic_data
    out: [full_user_data.json]
  assess_data:
    run: scripts/VIZ_prep_df.cwl
    label: 'Determine column types'
    in:
      data_frame: extend_user_formulas/full_user_data.json
    out: [continuous_cols, discrete_cols]
  fill_plot_options:
    run: ui-queries/any-viz.cwl
    label: 'Set plot options'
    doc: "Selections here pick the plot type and how it should be generated. For addtional details, see https://mountetna.github.io/vulcan.html#the-plot-configuration-interface which is clickably linked within this workflow's 'vignette'."
    in:
      data_frame: extend_user_formulas/full_user_data.json
      continuous_cols: assess_data/continuous_cols
      discrete_cols: assess_data/discrete_cols
    out: [plot_setup]
  make_plot:
    run: scripts/make_plot.cwl
    label: 'Create Plot'
    in:
      plot_setup: fill_plot_options/plot_setup
      data_frame: extend_user_formulas/full_user_data.json
    out: [plot.json, plot.png]
  show_plot:
    run: ui-outputs/plotly.cwl
    in:
      a: make_plot/plot.json
    out: []
    label: 'Display Plot'
