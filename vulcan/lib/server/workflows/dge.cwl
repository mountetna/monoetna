cwlVersion: v1.1
class: Workflow

inputs:
<<<<<<< HEAD
  1_Processed_data__link:
    type: str
    default: ''
    label: 'Provide a link?'
=======
  1_Cell_Filtering__by_covariate:
    type: int
    default: 1
    label: 'Filter cell by covariate'
    doc: 'TODO'
  2_DGE__dge:
    type: int
    default: 1
    label: 'TODO'
    doc: 'TODO'
  3_Select_output_types__outtypes:
    type: int
    default: 1
    label: 'TODO'
>>>>>>> parent of 38e0bf454 (Gvaihir - vulcan-dex (v. 0.0.1):)
    doc: 'TODO'

outputs:
  the_data:
    type: File
    outputSource: calc_umap/umap_anndata.h5ad

steps:
<<<<<<< HEAD
  get_data:
    run: scrpts/data_from_remote.cwl
    label: 'Provide a link to h5, cache h5'
    doc: 'TODO'
    in:
      url: 1_Processed_data__link
    out: [processed_data.h5ad,covariates.json]
  select_covariates:
    run: ui-queries/single_dropdown_multicheckbox.cwl
    label: 'Covariate selection'
    doc: 'TODO'
=======
  projectData:
    run: scripts/umap-projects.cwl
    label: 'Fetch project settings'
    in:
      a: 1_Cell_Filtering__by_covariate
      b: 2_DGE__dge
      c: 3_Select_output_types__outtypes
    out: [project_data]
  queryMagma:
    run: scripts/retrieve_selection_options.cwl
    label: 'Fetch selection options'
    in:
      project_data: projectData/project_data
    out: [selection_options]
  selectOnFeatures:
    run: ui-queries/multiple-multiselect-string-all.cwl
    label: 'Record Selection'
    doc: 'Selections here pick the subset of tube records to process and analyze. Select the values of the given features that you would like to target. The union of single-cell tube records that meet these criteria will be presented for confirmation, in the next step, based on the union of ALL feature selections here. Selections must be made for each in order to proceed, but if you want to just select tube records directly, pick the `All` option for all dropdowns here.'
>>>>>>> parent of 38e0bf454 (Gvaihir - vulcan-dex (v. 0.0.1):)
    in:
      opts: get_data/covariates.json
    out: [selected_options]
  perform_dge:
    run: scripts/dge.cwl
    label: 'Start actual DGE'
    doc: 'TODO'
    in:
      data: get_data/processed_data.h5ad
    out: [dge_data.h5ad]
  select_download_types:
    run: ui-queries/ceckboxes.cwl
    label: 'Select types of data to download'
    doc: 'TODO'
    in:
      data:
    out: [selected_outs.json]
  download_data:
    run: ui-outpits/link.ccwl
    label: 'Download links'
    doc: 'TODO'
    in:
      download_types: select_download_types/selected_outs.json
    out: []

