import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import StepName from './step_name';

import {
  validStep,
  hasUiInput,
  wrapEditableInputs,
  uiStepInputNames,
  uiStepType
} from '../../../utils/workflow';

export default function StepPending({step, stepIndex}) {
  const {workflow, pathIndex, session, status, setInputs} = useContext(
    VulcanContext
  );

  if (
    !validStep({workflow, pathIndex, stepIndex}) ||
    !session ||
    !status ||
    !step ||
    null === stepIndex ||
    !hasUiInput(step)
  )
    return null;

  function handleInputChange(inputName, value) {
    console.log('handling change', inputName, value);
    let userInputs = {...session.inputs};
    userInputs[inputName] = value;
    setInputs(userInputs);
  }

  // We need to map the user input step's output to
  //   a set of input items.
  let inputNames = uiStepInputNames(step);
  let mockStepInputs = inputNames.reduce((result, inputName) => {
    result[inputName] = {
      type: 'int', //uiStepType(step),
      label: step.label || step.name
    };
    return result;
  }, {});
  let components = wrapEditableInputs(mockStepInputs, handleInputChange);

  return (
    <div className='step-pending step'>
      <StepName
        step={step}
        status={status[pathIndex][stepIndex].status}
      ></StepName>
      <div className='step-inputs inputs-pane'>
        <div class='title'>User inputs</div>
        <div className='step-inputs-container items'>
          {components.map((comp) => {
            return comp;
          })}
        </div>
      </div>
    </div>
  );
}
