import React, {useContext} from 'react';
import {VulcanContext} from '../../../contexts/vulcan_context';
import StepIconName from './step_elements/step_icon_name';
import {statusOfStep} from '../../../selectors/workflow_selectors';
import {WorkspaceStep} from '../../../api_types';

export default function StepError({step}: {step: WorkspaceStep}) {
  let {state} = useContext(VulcanContext);
  const stepStatus = statusOfStep(step, state.status, state.workspace);

  if (!stepStatus) return null;

  let message = stepStatus.error || 'Something went wrong with this step.';

  return (
    <div className='step-error'>
      <StepIconName step={step} />
      <div className='text-wrapper'>
        <textarea readOnly rows={message.split('\n').length + 1}>
          {message}
        </textarea>
      </div>
    </div>
  );
}
