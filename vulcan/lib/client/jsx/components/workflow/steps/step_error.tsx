import React, {useContext} from 'react';
import {VulcanContext} from '../../../contexts/vulcan_context';
import StepName from './step_name';
import {statusOfStep} from '../../../selectors/workflow_selectors';
import {WorkflowStep} from '../../../api_types';

export default function StepError({step}: {step: WorkflowStep}) {
  let {state} = useContext(VulcanContext);
  const stepStatus = statusOfStep(step, state.status);

  if (!stepStatus) return null;

  let message = stepStatus.error || 'Something went wrong with this step.';

  return (
    <div className='step-error'>
      <StepName step={step} />
      <div className='text-wrapper'>
        <textarea readOnly rows={message.split('\n').length + 1}>
          {message}
        </textarea>
      </div>
    </div>
  );
}
