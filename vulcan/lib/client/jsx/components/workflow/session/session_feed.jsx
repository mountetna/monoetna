import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import {
  validPath,
  missingUiInputs,
  inputNamesToHashStub
} from '../../../utils/workflow';
import {
  completedUiStepsSelector,
  nextUiStepIndexSelector,
  errorStepsSelector,
  completedUiOutputsSelector
} from '../../../selectors/workflow';
import StepUserInput from '../steps/step_user_input';
import StepError from '../steps/step_error';
import StepOutput from '../steps/step_output';

export default function SessionFeed() {
  // Shows stream of Input, Output, Plots, etc.,
  //   as the session object updates.
  const context = useContext(VulcanContext);
  const {workflow, session, pathIndex, status, setInputs} = context;

  if (!workflow || !validPath({workflow, pathIndex}) || !session || !status)
    return null;

  let uiSteps = completedUiStepsSelector(context);
  let nextInputStepIndex = nextUiStepIndexSelector(context);
  let nextInputStep;

  // We inject a `null` input into the session,
  //   to indicate that we're waiting for a user input value
  //   to return to the server.
  if (-1 !== nextInputStepIndex) {
    nextInputStep = workflow.steps[pathIndex][nextInputStepIndex];

    let missingInputs = missingUiInputs(nextInputStep, session);

    if (missingInputs.length > 0) {
      setInputs(inputNamesToHashStub(missingInputs));
    }

    uiSteps.push({
      step: nextInputStep,
      index: nextInputStepIndex
    });
  }

  let errorSteps = errorStepsSelector(context);

  let outputs = completedUiOutputsSelector(context);

  return (
    <div className='session-feed'>
      {uiSteps.map((s, index) => (
        <StepUserInput
          key={index}
          step={s.step}
          stepIndex={s.index}
        ></StepUserInput>
      ))}
      {errorSteps.map((s, index) => (
        <StepError key={index} step={s.step} stepIndex={s.index}></StepError>
      ))}
      {outputs.map((s, index) => {
        <StepOutput key={index} step={s.step} stepIndex={s.index}></StepOutput>;
      })}
    </div>
  );
}
