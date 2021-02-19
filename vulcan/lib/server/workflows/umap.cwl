cwlVersion: v1.1
class: Workflow

inputs:
  someInt:
    type: int
    default: 200
    label: 'it is an int'
  someInt2:
    type: int
    default: 2
    label: 'another int'
  someIntWithoutDefault:
    type: int
    label: 'User supplies a number'

outputs:
  the_result:
    type: int
    outputSource: finalStep/sum

steps:
  firstAdd:
    run: scripts/add.cwl
    label: 'Add two numbers'
    in:
      a: someInt
      b: someInt2
    out: [sum]
  pickANum:
    run: ui-queries/multiselect.cwl
    in:
      a: firstAdd/sum
    out: [num]
  finalStep:
    run: scripts/add.cwl
    label: 'Add your number to the previous sum'
    in:
      a: firstAdd/sum
      b: pickANum/num
    out: [sum]
