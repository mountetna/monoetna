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

export default function SessionFeed() {
  // Shows stream of Input, Output, Plots, etc.,
  //   as the session object updates.
  const invoke = useActionInvoker();
  const context = useContext(VulcanContext);
  const {workflow, session, pathIndex, status} = context;

  if (!workflow || !validPath({workflow, pathIndex}) || !session || !status)
    return null;

  let stepsToRender = status[pathIndex]
    .map((step, index) => {
      if (STATUS.COMPLETE === step.status) {
        return {
          step: workflow.steps[pathIndex][index],
          index
        };
      }
    })
    .filter((s) => s);

  return (
    <div className='session-feed'>
      {stepsToRender.map((s) => (
        <StepComplete step={s.step} stepIndex={s.index}></StepComplete>
      ))}
    </div>
  );
}
