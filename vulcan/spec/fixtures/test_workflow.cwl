cwlVersion: v1.1
class: Workflow

inputs:
  someInt:
    type: int
    default: 200
    label: 'it is an int'
    doc: 'help tip'
  someIntWithoutDefault:
    type: int
    doc: 'another tip'

outputs:
  the_result:
    type: int
    outputSource: finalStep/sum
  thumbnail:
    type: File
    format: image/png
    outputSource: finalStep/thumb.png

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
    out: [sum, thumb.png]
  aPlot:
    run: ui-outputs/plotly.cwl
    in:
      a: finalStep/sum
      b: finalStep/thumb.png
    out: []
