import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import StepName from './step_name';

import {validStep} from '../../../utils/workflow';

export default function StepError({step, stepIndex}) {
  const {workflow, pathIndex, status} = useContext(VulcanContext);

  if (
    !validStep({workflow, pathIndex, stepIndex}) ||
    !status ||
    !step ||
    null === stepIndex
  )
    return null;

  let message =
    status[pathIndex][stepIndex].message ||
    'Something went wrong with this step.';

  return (
    <div className='step-error'>
      <StepName
        step={step}
        status={status[pathIndex][stepIndex].status}
      ></StepName>
      <div className='text-wrapper'>
        <textarea readonly rows={message.split('\n').length + 1}>
          {message}
        </textarea>
      </div>
    </div>
  );
}
