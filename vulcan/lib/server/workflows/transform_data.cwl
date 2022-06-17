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
  the_data:
    type: File
    outputSource: transform_data/formulaic_data

steps:
  get_data:
    run: scripts/VIZ_query_df.cwl
    label: 'Fetch Data'
    in:
      queryTerms: 1_Data_Source_magma_query__queryTerms
      user_columns: 1_Data_Source_magma_query__user_columns
      expand_matrices: 1_Data_Source_magma_query__expand_matrices
    out: [data_frame]
  transform_data:
    run: ui-queries/data-transformation.cwl
    label: 'Transform your data'
    doc: "Manipulate your data frame as needed. Right click for an interactions menu where you can add/remove columns. Start a cell with '=' to create functions in an Excel-like manner."
    in:
      data_frame: get_data/data_frame
    out: [formulaic_data]
  extend_user_formulas:
    run: scripts/calc_data_frame.cwl
    label: 'Extend formulas to all rows'
    in:
      original_data.json: get_data/data_frame
      user_data.json: transform_data/formulaic_data
    out: [full_user_data.json]
  show_data:
    run: ui-outputs/link.cwl
    in:
      a: extend_user_formulas/full_user_data.json
    out: []
    label: 'Download your final data frame'
  extra_show_data:
    run: ui-outputs/link.cwl
    in:
      a: extend_user_formulas/full_user_data.json
    out: []
    label: 'Download your final data frame again!'