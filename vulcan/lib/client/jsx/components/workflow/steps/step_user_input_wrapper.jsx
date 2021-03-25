import React, {useContext, useState, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import {STATUS} from '../../../models/steps';
import StepName from './step_name';
import StepUserInput from './step_user_input';
import StepUserInputDrawer from './step_user_input_drawer';

import {
  validStep,
  hasUiInput,
  removeDependentInputs,
  dependentInputsData,
  isMissingDrawerInputs,
  isMissingStandardInputs
} from '../../../utils/workflow';

export default function StepUserInputWrapper({step, stepIndex}) {
  const [open, setOpen] = useState(true);
  const {workflow, pathIndex, session, status, setInputs, setData} = useContext(
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
    // Also have to set the downstream steps' inputs
    //   data to null, so the items are removed from
    //   the InputFeed.
    dependentInputsData({
      workflow,
      status,
      inputName: inputName.split('/')[0]
    }).forEach((action) => {
      setData(action.url, action.data);
    });

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

    if (step.out && isMissingStandardInputs(step, session)) {
      setOpen(true);
    } else if (!step.out && isMissingDrawerInputs(step, session)) {
      setOpen(true);
    } else if (STATUS.COMPLETE === stepStatus) {
      setOpen(false);
    }
  }, [status, session]);

  let Component = step.isGroup ? StepUserInputDrawer : StepUserInput;

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
        <Component
          key={`${stepIndex}-${status[pathIndex][stepIndex].status}`}
          step={step}
          handleInputChange={handleInputChange}
        ></Component>
      </div>
    </div>
  );
}
