cwlVersion: v1.1
class: Workflow

inputs: []

outputs:
  the_csv:
    type: File
    outputSource: DGEcalc/diffexp.csv

steps:
  mockSetup:
    run: scripts/dge_mock_data.cwl
    label: 'Fetch data'
    in: []
    out: [scdata.h5ad]
    
  parse_metadata_options_full:
    run: scripts/dge_summarize_obs.cwl
    label: 'Determine Subsetting Options'
    in:
      scdata.h5ad: mockSetup/scdata.h5ad
    out: [discrete_metas, continuous_metas]
    
  selectSubsetMeta:
    run: ui-queries/select-autocomplete.cwl
    label: 'Select metadata to subset on'
    doc: 'Selections here pick the attribute with which to group cells for subsetting which cells DGE will be run on. Actual group-names (values for this data) will be presented in the next step.'
    in:
      a: parse_metadata_options_full/discrete_metas
    out: [subset_meta]
    
  prepSubsetVals:
    run: scripts/dge_retrieve_meta_vals.cwl
    label: 'Parse metadata selection'
    in:
      scdata.h5ad: mockSetup/scdata.h5ad
      meta: selectSubsetMeta/subset_meta
    out: [opts]
    
  selectSubsetVals:
    run: ui-queries/multiselect-string.cwl
    label: 'Select Subset Targets'
    doc: 'Selections here pick the subset of attributes base cell subsetting on. Actual values for this data will be presented in the next step. If you DO NOT WISH TO SUBSET to certain cells, select any one of these, then "all" in the next step.'
    in:
      a: prepSubsetVals/opts
    out: [selected_vals]
    
  subsetData:
    run: scripts/dge_subset.cwl
    label: 'Subset cells'
    in:
      selected_vals: selectSubsetVals/selected_vals
      meta: selectSubsetMeta/subset_meta
      scdata.h5ad: mockSetup/scdata.h5ad
    out: [scdata.h5ad, pauser, methods]
    
  selectDGEMeta:
    run: ui-queries/select-autocomplete.cwl
    label: 'Select metadata to run DGE on'
    doc: 'Selections here pick the attribute with which to group cells for differential expression analysis. Actual group-names (values for this data) will be presented in the next step.'
    in:
      a: parse_metadata_options_full/discrete_metas
      b: subsetData/pauser
    out: [dge_meta]
    
  prepDGEVals:
    run: scripts/dge_retrieve_meta_vals.cwl
    label: 'Parse metadata selection'
    in:
      scdata.h5ad: subsetData/scdata.h5ad
      meta: selectDGEMeta/dge_meta
    out: [opts]
    
  selectDGEVals_1:
    run: ui-queries/multiselect-string.cwl
    label: 'Select DGE Target-1'
    in:
      a: prepDGEVals/opts
    out: [selected_vals]
    
  selectDGEVals_2:
    run: ui-queries/multiselect-string.cwl
    label: 'Select DGE Target-2'
    in:
      a: prepDGEVals/opts
    out: [selected_vals]
    
  DGEmethod:
    run: ui-queries/select-autocomplete.cwl
    doc: "Your selection here picks the 'method' of the 'scanpy.tl.rank_genes_groups' differential expression function."
    label: 'Select p-value Calculation Method'
    in:
      a: subsetData/methods
    out: [dge_method]
  
  DGEcalc:
    run: scripts/dge_calc.cwl
    label: 'Calculate DGE'
    in:
      scdata.h5ad: subsetData/scdata.h5ad
      target_meta: selectDGEMeta/dge_meta
      target_1: selectDGEVals_1/selected_vals
      target_2: selectDGEVals_2/selected_vals
      dge_method: DGEmethod/dge_method
    out: [diffexp.csv]
    
  downloadDEData:
    run: ui-outputs/link.cwl
    in:
      a: DGEcalc/diffexp.csv
    out: []
    label: 'Download DiffExp as csv'
