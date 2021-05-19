cwlVersion: v1.1
class: Workflow

inputs:
  1_DGE_Calculation_Options__dge_method:
    type: string
    default: 'wilcoxon'
    label: 'testing method'
    doc: 'A string indicating what scanpy method option to use for calculating differential expression. Options are: "logreg", "t-test", "wilcoxon", "t-test_overestim_var". See documentation for "scanpy.tl.rank_genes_groups" for further details.'
  HOLDER__umap_session_info:
    type: string
    label: 'FILLER'
    default: 'FILLER'

outputs:
  the_data:
    type: File
    outputSource: DGEcalc/scdata
  the_csv:
    type: File
    outputSource: DGEcalc/diffexp.csv

steps:
  projectData:
    run: scripts/umap-projects.cwl
    label: 'Fetch project settings'
    in: []
    out: [project_data]
  queryMagma:
    run: scripts/retrieve_aattribute_options.cwl
    label: 'Fetch data options'
    in:
      project_data: projectData/project_data
      umap_session_info: HOLDER_umap_session_info
    out: [selection_atts, scData]
  selectSubsetAttributes:
    run: ui-queries/multiselect-string-all.cwl
    label: 'Select Subset Targets-1'
    doc: 'Selections here pick the subset of attributes base cell subsetting on. Actual values for this data will be presented in the next step. If you DO NOT WISH TO SUBSET to certain cells, select any one of these, then "all" in the next step.'
    in:
      a: queryMagma/selection_atts
    out: [selected_atts]
  parseSelectionAttributes:
    run: scripts/parse_attribute_selections.cwl
    label: 'Parse Subset Targets-1'
    in:
      project_data: projectData/project_data
      selected_atts: selectSubsetAttributes/selected_atts
      scdata: queryMagma/scData
    out: [selection_values]
  selectSubsetValues:
    run: ui-queries/multiple-multiselect-string-all.cwl
    label: 'Select Subseting Targets-2'
    doc: 'Selections here pick the subset of cells to move forward with. Select all values of the selected attributes that you would like to target. The cells that meet these criteria, based on the union of ALL feature selections here, will be retained for differential expression calculation steps. Selections must be made for each in order to proceed. If you DO NOT WISH TO SUBSET to certain cells, select "all".'
    in:
      a: parseSelectionAttributes/selection_values
    out: [selected_values]
  subsetData:
    run: scripts/parse_value_selections.cwl
    label: 'Subset cells'
    in:
      scdata: queryMagma/scdata
      selected_values: selectSubsetValues/selected_values
    out: [scdata]
  selectDGEAttributes:
    run: ui-queries/multiselect-string-all.cwl
    label: 'Select DGE Target-1'
    doc: 'Selections here pick the attribute with which to group cells for differential expression analysis. Actual group-names (values for this data) will be presented in the next step.'
    in:
      a: queryMagma/selection_atts
    out: [selected_att]
  parseDGEAttributes:
    run: scripts/parse_attribute_selections.cwl
    label: 'Parse DGE Target-1'
    in:
      selected_att: selectDGEAttributes/selected_att
      scdata: subsetData/scData
    out: [selection_values]
  selectDGEGroups:
    run: ui-queries/multiselect-string-all.cwl
    label: 'Select DGE Groups'
    doc: 'Selections here pick the groups to compare in DGE analysis. In an initial iteration of this workflow, cells of each individual selection here will be compared to all other cells.'
    in:
      a: parseDGEAttributes/selection_values
    out: [dge_groups]
  DGEcalc:
    run: scripts/dge_calc.cwl
    label: 'Calculate DGE'
    in:
      scdata: subsetData/scData
      target_att: selectDGEAttributes/selected_att
      target_groups: selectDGEGroups/dge_groups
      dge_method: 1_DGE_Calculation_Options__dge_method
    out: [scdata, diffexp.csv]
  downloadDEData:
    run: ui-outputs/link.cwl
    in:
      a: DGEcalc/diffexp.csv
    out: []
    label: 'Download DiffExp as csv'

