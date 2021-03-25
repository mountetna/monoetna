import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import {
  validPath,
  missingUiInputs,
  inputNamesToHashStub,
  groupUiSteps
} from '../../../utils/workflow';
import {
  completedUiStepsSelector,
  nextUiStepsSelector,
  errorStepsSelector
} from '../../../selectors/workflow';

import PrimaryInputs from './primary_inputs';
import StepUserInputWrapper from '../steps/step_user_input_wrapper';
import StepError from '../steps/step_error';

export default function InputFeed() {
  // Shows stream of Inputs,
  //   as the session object updates.
  const context = useContext(VulcanContext);
  const {workflow, session, pathIndex, status, setInputs} = context;

  if (!workflow || !validPath({workflow, pathIndex}) || !session || !status)
    return null;

  let uiSteps = completedUiStepsSelector(context);
  let nextUiSteps = nextUiStepsSelector(context);

  // We inject a `null` input into the session,
  //   to indicate that we're waiting for a user input value
  //   to return to the server.
  if (nextUiSteps.length > 0) {
    nextUiSteps.forEach((nextInputStep) => {
      let missingInputs = missingUiInputs(nextInputStep, session);

      if (missingInputs.length > 0) {
        console.log('with some missing inputs', missingInputs);
        // Make sure to copy over the current inputs, otherwise
        //   they'll get wiped out in the reducer.
        setInputs({
          ...session.inputs,
          ...inputNamesToHashStub(missingInputs)
        });
      }

      uiSteps.push({
        step: nextInputStep,
        index: nextInputStep.index
      });
    });
  }

  let errorSteps = errorStepsSelector(context);

  uiSteps = groupUiSteps(uiSteps);

  return (
    <div className='session-input-feed'>
      <PrimaryInputs></PrimaryInputs>
      {uiSteps.map((s, index) => (
        <StepUserInputWrapper
          key={index}
          step={s.step}
          stepIndex={s.index}
        ></StepUserInputWrapper>
      ))}
      {errorSteps.map((s, index) => (
        <StepError key={index} step={s.step} stepIndex={s.index}></StepError>
      ))}
    </div>
  );
}
