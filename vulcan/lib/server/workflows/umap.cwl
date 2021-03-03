cwlVersion: v1.1
class: Workflow

inputs:
  includeParamA:
    type: boolean
    default: true
    label: 'regress by param A?'
  includeParamB:
    type: boolean
    default: true
    label: 'regress by param B?'
  includeParamC:
    type: boolean
    default: true
    label: 'regress by param C?'
  maxPcs:
    type: int
    default: 15
    label: 'Maximum number of PCs'

outputs:
  the_data:
    type: File
    outputSource: doUmapStuff/umap

steps:
  queryMagma:
    run: scripts/fake_query.cwl
    label: 'Fetch pool record names'
    in:
      b: includeParamA
      c: includeParamB
      d: includeParamC
      e: maxPcs
    out: [names]
  pickPools:
    run: ui-queries/multiselect-string.cwl
    label: 'Select pool records'
    in:
      a: queryMagma/names
    out: [names]
  doUmapStuff:
    run: scripts/fake_large_umap_calculations.cwl
    label: 'Do the UMAP stuff'
    in:
      a: pickPools/names
      b: includeParamA
      c: includeParamB
      d: includeParamC
      e: maxPcs
    out: [umap]
  showMe:
    run: ui-outputs/plotly.cwl
    in:
      a: doUmapStuff/umap
    out: []
    label: 'UMAP'
  convertToConsignment:
    run: scripts/convert_to_consignment.cwl
    label: 'Convert data to TSV'
    in:
      a: doUmapStuff/umap
    out: [output]
  downloadRawData:
    run: ui-outputs/consignment.cwl
    in:
      a: convertToConsignment/output
    out: []
    label: 'Download Raw Data as TSV'