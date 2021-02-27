import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import {validStep, wrapPaneItem, stringify} from '../../../utils/workflow';

import StepName from './step_name';

// Small, non-editable version of a
//   step, used for informational / non-interactive purposes.
export default function StepViewCard({step, stepIndex}) {
  const {workflow, pathIndex, session, status} = useContext(VulcanContext);

  if (
    !validStep({workflow, pathIndex, stepIndex}) ||
    !session ||
    !status ||
    !step ||
    null === stepIndex
  )
    return null;

  // Find the step's inputs + values from the Primary Inputs
  let inputValues = step.in.map((input) => {
    let value;
    let name = input.id;
    if ('primary_inputs' === input.source[0]) {
      let primaryInputName = input.source[1];
      value = session.inputs[primaryInputName];
    } else {
      let outputStepName = input.source[0];
      let outputStepIndex = workflow.steps[pathIndex].findIndex(
        (s) => outputStepName === s.name
      );
      let outputVariableName = input.source[1];
      // Sometimes data won't be available yet, so we
      //   have to punt and wait for State to update.
      if (
        status[pathIndex][outputStepIndex] &&
        status[pathIndex][outputStepIndex].data
      ) {
        value = stringify(
          status[pathIndex][outputStepIndex].data[outputVariableName]
        );
      }
    }

    return {
      name,
      value
    };
  });

  // For the output values, we'll have to figure out the type
  //   of thing to render. For now, just show the raw data.
  let outputValue;
  let outputName = step.out[0];

  if (
    status[pathIndex][stepIndex].data &&
    status[pathIndex][stepIndex].data[outputName]
  ) {
    outputValue = stringify(status[pathIndex][stepIndex].data[outputName]);
  }

  return (
    <div className='step-view-card step'>
      <StepName
        step={step}
        status={status[pathIndex][stepIndex].status}
      ></StepName>
      <div className='step-inputs inputs-pane'>
        <div class='title'>Step inputs</div>
        <div className='step-inputs-container items'>
          {inputValues.map((input, index) => {
            return wrapPaneItem(input, index);
          })}
        </div>
      </div>
      <div className='step-outputs outputs-pane'>
        <div class='title'>Step outputs</div>
        <div className='step-outputs-container items'>
          {wrapPaneItem({
            name: outputName,
            value: outputValue
          })}
        </div>
      </div>
    </div>
  );
}
