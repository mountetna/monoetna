cwlVersion: v1.1
class: Workflow

inputs:
  1_Data_Input__string:
    type: string
    label: "String"
  1_Data_Input__number:
    type: float
    label: "Number"
  1_Data_Input__integer:
    type: int
    label: "Integer"
  1_Data_Input__boolean:
    type: boolean
    label: "Boolean"
    default: true

outputs:
  the_data:
    type: File
    outputSource: get_data/str

steps:
  get_data:
    run: scripts/test_R_dataflow.cwl
    label: 'Read and Write in R'
    in:
      string: 1_Data_Input__string
      number: 1_Data_Input__number
      integer: 1_Data_Input__integer
      boolean: 1_Data_Input__boolean
    out: [str,num,int,bool]
  download_string:
    run: ui-outputs/link.cwl
    in:
      a: filter_DGE/filtered_diffexp.csv
    out: []
    label: 'Download String output'
  download_number:
    run: ui-outputs/link.cwl
    in:
      a: get_data/str
    out: []
    label: 'Download Number output'
  download_boolean:
    run: ui-outputs/link.cwl
    in:
      a: get_data/int
    out: []
    label: 'Download Integer output'
  download_boolean:
    run: ui-outputs/link.cwl
    in:
      a: get_data/bool
    out: []
    label: 'Download Boolean output'