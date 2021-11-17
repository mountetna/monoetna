import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {WorkflowStep} from '../../../api_types';
import {labelOfStepOrGroupedStep} from '../../../selectors/workflow_selectors';
import {STATUS} from '../../../api_types';
import {WorkflowStepGroup} from "../user_interactions/inputs/input_types";

import StepIcon from './step_icon';

const StepName = ({
  step
}: {
  step: WorkflowStep | WorkflowStepGroup
}) => {
  let {state} = useContext(VulcanContext);
  const label = labelOfStepOrGroupedStep(step);

  return (
    <div className='step-name'>
      <StepIcon step={step}/>
      <div className='step-button'>{label}</div>
    </div>
  );
};

export default StepName;
