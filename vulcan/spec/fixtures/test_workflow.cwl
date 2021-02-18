cwlVersion: v1.1
class: Workflow

inputs:
  someInt:
    type: int
    default: 200
    label: 'it is an int'
  someIntWithoutDefault:
    type: int

outputs:
  the_result:
    type: int
    outputSource: finalStep/sum

steps:
  firstAdd:
    run: scripts/add.cwl
    in:
      a: someInt
      b: someIntWithoutDefault
    out: [sum]
  pickANum:
    run: ui-queries/pick-a-number.cwl
    in:
      num: firstAdd/sum
    out: [num]
  finalStep:
    run: scripts/add.cwl
    in:
      a: firstAdd/sum
      b: pickANum/num
    out: [sum]
