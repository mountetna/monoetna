cwlVersion: v1.1
class: Workflow

inputs:
  1_Processed_data__link:
    type: string
    default: ''
    label: 'Provide a link?'
    doc: 'TODO'

outputs:
  the_data:
    type: File
    outputSource: dge/processed_anndata.h5ad

steps:
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
    in:
      opts: get_data/covariates.json
    out: [selected_options.json]
  perform_dge:
    run: scripts/dge.cwl
    label: 'Start actual DGE'
    doc: 'TODO'
    in:
      selected_covariates: select_covariates/selected_options.json
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



