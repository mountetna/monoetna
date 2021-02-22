import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import {validPath} from '../../../selectors/workflow_selector';

import StepName from './step_name';

export default function StepComplete({step, stepIndex}) {
  const {workflow, pathIndex, session, status} = useContext(VulcanContext);

  if (
    !validPath({workflow, pathIndex}) ||
    !session ||
    !status ||
    !step ||
    null === stepIndex
  )
    return null;

  function wrapPaneItem(item) {
    return (
      <div className='view_item'>
        <div className='item_name'>{item.name}</div>
        <div className='item_view'>{item.value}</div>
      </div>
    );
  }

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
      value = status[(pathIndex, outputStepIndex)].data[outputVariableName];
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
    outputValue = status[pathIndex][stepIndex].data[outputName];
  }

  return (
    <div className='step-complete step'>
      <StepName
        step={step}
        status={status[pathIndex][stepIndex].status}
      ></StepName>
      <div className='step-inputs inputs-pane'>
        <div class='title'>Step inputs</div>
        <div className='step-inputs-container items'>
          {inputValues.map((input) => {
            return wrapPaneItem(input);
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
