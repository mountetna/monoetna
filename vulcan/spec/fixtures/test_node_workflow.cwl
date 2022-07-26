cwlVersion: v1.1
class: Workflow

inputs:
  someInt:
    type: int
    default: 200
    label: 'it is an int'
    doc: 'help tip'

outputs:
  the_result:
    type: string
    outputSource: makeDF/df.json

steps:
  makeDF:
    run: scripts/mock_data_frame.cwl
    in:
      a: someInt
    out: [df.json]

