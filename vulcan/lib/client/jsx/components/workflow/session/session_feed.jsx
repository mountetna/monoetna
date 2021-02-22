import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import Icon from 'etna-js/components/icon';

import {submit} from '../../../api/vulcan';
import {VulcanContext} from '../../../contexts/vulcan';
import {
  allInputsDefined,
  defaultInputValues,
  validPath
} from '../../../selectors/workflow_selector';
import {STATUS} from '../../../models/steps';
import StepComplete from '../steps/step_complete';
import StepPending from '../steps/step_pending';

export default function SessionFeed() {
  // Shows stream of Input, Output, Plots, etc.,
  //   as the session object updates.
  const invoke = useActionInvoker();
  const context = useContext(VulcanContext);
  const {workflow, session, pathIndex, status, setInputs} = context;

  if (!workflow || !validPath({workflow, pathIndex}) || !session || !status)
    return null;

  let completedSteps = status[pathIndex]
    .map((step, index) => {
      if (STATUS.COMPLETE === step.status) {
        return {
          step: workflow.steps[pathIndex][index],
          index
        };
      }
    })
    .filter((s) => s);

  let nextInputStepIndex = status[pathIndex].findIndex(
    (s) => STATUS.PENDING === s.status
  );

  // We need to inject a `null` input into the session,
  //   to indicate that we're waiting for a user input value
  //   to return to the server.
  if (-1 !== nextInputStepIndex) {
    let nextInputStep = workflow.steps[pathIndex][nextInputStepIndex];
    let missingInputs = nextInputStep.in.filter((input) => {
      return !Object.keys(session.inputs).includes(input.id);
    });
    setInputs(
      missingInputs.reduce((result, input) => {
        if (!result.hasOwnProperty(inptu.id)) {
          result[input.id] = null;
        }
      }, {})
    );
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
