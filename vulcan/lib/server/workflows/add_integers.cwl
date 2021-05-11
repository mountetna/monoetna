cwlVersion: v1.1
class: Workflow

inputs:
  anotherInt:
    type: int
    default: 500
    label: 'A'
    doc: 'Just put in an integer'
  someInt:
    type: int
    default: 200
    label: 'B'

outputs:
  the_result:
    type: int
    outputSource: finalStep/sum

steps:
  firstAdd:
    run: scripts/add.cwl
    in:
      a: someInt
      b: anotherInt
    out: [sum]
    label: 'Add two integers'
  pickANum:
    run: ui-queries/int.cwl
    in:
      num: firstAdd/sum
    out: [num]
    label: 'Provide an integer'
  finalStep:
    run: scripts/add.cwl
    in:
      a: firstAdd/sum
      b: pickANum/num
    out: [sum]
    label: 'Final sum'
  showFinalSum:
    run: ui-outputs/raw.cwl
    in:
      a: finalStep/sum
    out: []
    label: 'View your sum'
