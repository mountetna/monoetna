import React, {useContext, useState, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import {STATUS} from '../../../models/steps';
import UserInput from '../user_interactions/inputs/user_input';
import StepName from './step_name';

import {
  validStep,
  hasUiInput,
  uiStepInputNames,
  uiStepType,
  uiStepOptions,
  removeDependentInputs
} from '../../../utils/workflow';

export default function StepUserInput({step, stepIndex}) {
  const [open, setOpen] = useState(true);
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

    // Changing a user input should remove any
    //   subsequent inputs from userInputs,
    //   because we don't know how the selection
    //   changes downstream calculations.
    // If we leave them in, the workflow will
    //   re-run with potentially invalid inputs
    //   to subsequent steps.
    setInputs(
      removeDependentInputs({
        userInputs,
        inputName: inputName.split('/')[0],
        workflow
      })
    );
  }

  function toggleInputs() {
    setOpen(!open);
  }

  useEffect(() => {
    let stepStatus = status[pathIndex][stepIndex].status;
    if (STATUS.COMPLETE === stepStatus) {
      setOpen(false);
    } else if (STATUS.PENDING === stepStatus) {
      setOpen(true);
    }
  }, [status]);

  // We need to map the user input step's output to
  //   a set of input items.
  let inputNames = uiStepInputNames(step);
  let stepInputs = inputNames.reduce((result, inputName) => {
    result.push({
      type: uiStepType(step),
      label: step.label || step.name,
      default: session.inputs[inputName] || null,
      options: uiStepOptions({step, pathIndex, status}),
      name: inputName
    });
    return result;
  }, []);

  return (
    <div className='step-user-input step'>
      <div onClick={toggleInputs}>
        <StepName
          step={step}
          status={status[pathIndex][stepIndex].status}
        ></StepName>
      </div>
      <div
        className={`step-user-input-inputs sliding-panel vertical ${
          open ? 'open' : 'closed'
        }`}
      >
        {stepInputs.map((input, index) => {
          return (
            <UserInput
              input={input}
              hideLabel={true}
              onChange={handleInputChange}
              key={index}
            ></UserInput>
          );
        })}
      </div>
    </div>
  );
}
