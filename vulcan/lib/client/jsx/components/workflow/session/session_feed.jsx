import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import {validPath, hasUiInput, uiStepInputNames} from '../../../utils/workflow';
import {STATUS} from '../../../models/steps';
import StepUserInput from '../steps/step_user_input';
import StepError from '../steps/step_error';

export default function SessionFeed() {
  // Shows stream of Input, Output, Plots, etc.,
  //   as the session object updates.
  const context = useContext(VulcanContext);
  const {workflow, session, pathIndex, status, setInputs} = context;

  if (!workflow || !validPath({workflow, pathIndex}) || !session || !status)
    return null;

  let uiSteps = status[pathIndex]
    .map((step, index) => {
      let workflowStep = workflow.steps[pathIndex][index];
      if (STATUS.COMPLETE === step.status && hasUiInput(workflowStep)) {
        return {
          step: workflowStep,
          index
        };
      }
    })
    .filter((s) => s);

  let nextInputStepIndex = status[pathIndex].findIndex((s, index) => {
    let workflowStep = workflow.steps[pathIndex][index];
    return STATUS.PENDING === s.status && hasUiInput(workflowStep);
  });
  let nextInputStep;

  // We inject a `null` input into the session,
  //   to indicate that we're waiting for a user input value
  //   to return to the server.
  if (-1 !== nextInputStepIndex) {
    nextInputStep = workflow.steps[pathIndex][nextInputStepIndex];

    let missingInputs = uiStepInputNames(nextInputStep).filter(
      (outputName) => !Object.keys(session.inputs).includes(outputName)
    );

    if (missingInputs.length > 0) {
      let missingInputsHash = missingInputs.reduce((result, input) => {
        if (!result.hasOwnProperty(input)) {
          result[input] = null;
        }

        return result;
      }, {});

      setInputs(missingInputsHash);
    }

    uiSteps.push({
      step: nextInputStep,
      index: nextInputStepIndex
    });
  }

  let errorSteps = status[pathIndex]
    .map((step, index) => {
      let workflowStep = workflow.steps[pathIndex][index];
      if (STATUS.ERROR === step.status) {
        return {
          step: workflowStep,
          index
        };
      }
    })
    .filter((s) => s);

  return (
    <div className='session-feed'>
      {uiSteps.map((s) => (
        <StepUserInput step={s.step} stepIndex={s.index}></StepUserInput>
      ))}
      {errorSteps.map((s) => (
        <StepError step={s.step} stepIndex={s.index}></StepError>
      ))}
    </div>
  );
}
