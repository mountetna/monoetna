import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import StepName from './step_name';

import {
  validStep,
  hasUiInput,
  wrapEditableInputs,
  uiStepInputNames,
  uiStepType,
  uiStepOptions
} from '../../../utils/workflow';

export default function StepUserInput({step, stepIndex}) {
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
    let userInputs = {...session.inputs};
    userInputs[inputName] = value;
    setInputs(userInputs);
  }

  // We need to map the user input step's output to
  //   a set of input items.
  let inputNames = uiStepInputNames(step);
  let mockStepInputs = inputNames.reduce((result, inputName) => {
    result[inputName] = {
      type: uiStepType(step),
      label: step.label || step.name,
      default: session.inputs[inputName] || null,
      options: uiStepOptions({step, pathIndex, status})
    };
    return result;
  }, {});
  let components = wrapEditableInputs(mockStepInputs, handleInputChange);

  return (
    <div className='step-user-input step'>
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
