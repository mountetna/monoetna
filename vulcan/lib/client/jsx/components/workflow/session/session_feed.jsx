import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';
import {validPath, hasUiInput, uiStepInputNames} from '../../../utils/workflow';
import {STATUS} from '../../../models/steps';
import StepComplete from '../steps/step_complete';
import StepPending from '../steps/step_pending';

export default function SessionFeed() {
  // Shows stream of Input, Output, Plots, etc.,
  //   as the session object updates.
  const context = useContext(VulcanContext);
  const {workflow, session, pathIndex, status, setInputs} = context;

  if (!workflow || !validPath({workflow, pathIndex}) || !session || !status)
    return null;

  let completedSteps = status[pathIndex]
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

  let nextInputStepIndex = status[pathIndex].findIndex(
    (s) => STATUS.PENDING === s.status
  );
  let nextInputStep;

  // TODO: This needs to be seriously refactored!
  // We inject a `null` input into the session,
  //   to indicate that we're waiting for a user input value
  //   to return to the server.
  if (-1 !== nextInputStepIndex) {
    nextInputStep = workflow.steps[pathIndex][nextInputStepIndex];

    if (hasUiInput(nextInputStep)) {
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
    }
  }

  return (
    <div className='session-feed'>
      {completedSteps.map((s) => (
        <StepComplete step={s.step} stepIndex={s.index}></StepComplete>
      ))}
      {nextInputStep ? (
        <StepPending
          step={nextInputStep}
          stepIndex={nextInputStepIndex}
        ></StepPending>
      ) : null}
    </div>
  );
}
