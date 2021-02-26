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
  includeParamA:
    type: boolean
    default: true
    label: 'regress by param A?'
  aFloat:
    type: float

outputs:
  the_result:
    type: int
    outputSource: finalStep/sum
    format: text

steps:
  firstAdd:
    run: scripts/add.cwl
    label: 'Add two numbers'
    in:
      a: someInt
      b: someInt2
    out: [sum]
  pickANum:
    run: ui-queries/multiselect-string.cwl
    in:
      a: firstAdd/sum
    out: [num]
  finalStep:
    run: scripts/add_from_list.cwl
    label: 'Add your numbers to the previous sum'
    in:
      a: firstAdd/sum
      b: pickANum/num
    out: [sum]