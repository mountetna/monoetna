cwlVersion: v1.1
class: Workflow

inputs:
  a:
    type: int
  b:
    type: int

outputs:
  the_result:
    type: int
    outputSource: finalStep/sum

steps:
  firstAdd:
    run: scripts/add.cwl
    in:
      a: a
      b: b
    out: [sum]
  pickANum:
    run: ui-queries/pick-a-number.cwl
    in: []
    out: [num]
  otherPickANumber:
    run: ui-queries/pick-a-number.cwl
    in: []
    out: [num]
  twiceOtherPickANumber:
    run: scripts/add.cwl
    in:
      a: otherPickANumber/num
      b: otherPickANumber/num
    out: [sum]
  twiceOtherPickANumber2:
    run: scripts/add.cwl
    in:
      a: otherPickANumber/num
      b: otherPickANumber/num
    out: [sum]
  thriceOtherPickANumber:
    run: scripts/add.cwl
    in:
      a: otherPickANumber/num
      b: twiceOtherPickANumber/sum
    out: [sum]
  ficeOtherPickANumber:
    run: scripts/add.cwl
    in:
      a: thriceOtherPickANumber/sum
      b: twiceOtherPickANumber2/sum
    out: [sum]
  finalStep:
    run: scripts/add.cwl
    in:
      a: pickANum/num
      b: ficeOtherPickANumber/sum
    out: [sum]
  aPlot:
    run: ui-outputs/plotter.cwl
    in:
      a: finalStep/sum
      b: firstAdd/sum
    out: []
