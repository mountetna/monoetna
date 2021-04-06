import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import Icon from 'etna-js/components/icon';
import AnimatedClock from './animated_clock';
import {WorkflowStep} from "../../../api_types";
import {statusOfStep, statusStringOfStepOrGroupedStep} from "../../../selectors/workflow_selectors";
import {STATUS} from "../../../api_types";
import {GroupedInputStep} from "../user_interactions/inputs/input_types";

const icons: {[k: string]: {icon: string, className: string}} = {
  [STATUS.COMPLETE]: { icon: 'check', className: 'light green' },
  [STATUS.RUNNING]: {icon: 'clock', className: 'light'},
  [STATUS.PENDING]: {icon: 'clock', className: 'light'},
  [STATUS.ERROR]: {icon: 'times-circle', className: 'light red'},
};

export default function StepName({step}: {step: WorkflowStep | GroupedInputStep}) {
  let {state, statusIsFresh} = useContext(VulcanContext);
  const {workflow, status} = state;
  if (!workflow) return null;
  const statusStr = statusStringOfStepOrGroupedStep(step, workflow, status);
  let icon = icons[statusStr] || icons[STATUS.PENDING];

  let className = `step-status-icon ${icon.className}`;
  let IconComponent = <Icon title={ step.label || step.name } className={className} icon={icon.icon}/>;

  if (STATUS.RUNNING === statusStr && !statusIsFresh) {
    IconComponent = <AnimatedClock />;
  }

  return (
    <div className='step-name'>
      <div className='step-status-icon-wrapper'>{IconComponent}</div>
      <div className='step-button'>{step.label || step.name}</div>
    </div>
  );
}
