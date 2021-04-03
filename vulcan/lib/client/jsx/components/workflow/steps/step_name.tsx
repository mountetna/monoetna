import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import Icon from 'etna-js/components/icon';
import AnimatedClock from './animated_clock';
import {WorkflowStep} from "../../../api_types";
import {statusOfStep} from "../../../selectors/workflow_selectors";
import {STATUS} from "../../../api_types";

const icons: {[k: string]: {icon: string, className: string}} = {
  [STATUS.COMPLETE]: { icon: 'check', className: 'light green' },
  [STATUS.RUNNING]: {icon: 'clock', className: 'light'},
  [STATUS.PENDING]: {icon: 'clock', className: 'light'},
  [STATUS.ERROR]: {icon: 'times-circle', className: 'light red'},
};

export default function StepName({step}: {step: WorkflowStep}) {
  let {state, isLoading} = useContext(VulcanContext);
  const stepStatus = statusOfStep(step, state.status);
  const statusStr = stepStatus ? stepStatus.status : STATUS.PENDING;
  let icon = icons[statusStr] || icons[STATUS.PENDING];

  let className = `step-status-icon ${icon.className}`;
  let IconComponent = <Icon className={className} icon={icon.icon}/>;

  if (STATUS.RUNNING === statusStr && isLoading) {
    IconComponent = <AnimatedClock />;
  }

  return (
    <div className='step-name'>
      <div className='step-status-icon-wrapper'>{IconComponent}</div>
      <div className='step-button'>{step.label || step.name}</div>
    </div>
  );
}
