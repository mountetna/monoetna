import React from 'react';

import {WorkspaceStep} from '../../../api_types';
import {labelOfStepOrGroupedStep} from '../../../selectors/workflow_selectors';
import {WorkspaceStepGroup} from '../user_interactions/inputs/input_types';

import StepIcon from './step_icon';

const StepName = ({
  step
}: {
  step: WorkspaceStep | WorkspaceStepGroup
}) => {
  const label = labelOfStepOrGroupedStep(step);

  return (
    <div className='step-name'>
      <StepIcon step={step}/>
      <div className='step-button'>{label}</div>
    </div>
  );
};

export default StepName;
