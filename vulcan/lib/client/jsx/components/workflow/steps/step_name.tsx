import React, {useContext} from 'react';

import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ExpandLessIcon from '@material-ui/icons/ExpandLess';

import {VulcanContext} from '../../../contexts/vulcan_context';
import Icon from 'etna-js/components/icon';
import AnimatedClock from './animated_clock';
import {WorkflowStep} from '../../../api_types';
import {statusStringOfStepOrGroupedStep} from '../../../selectors/workflow_selectors';
import {STATUS} from '../../../api_types';
import {GroupedInputStep} from '../user_interactions/inputs/input_types';

const icons: {[k: string]: {icon: string; className: string}} = {
  [STATUS.COMPLETE]: {icon: 'check', className: 'light green'},
  [STATUS.RUNNING]: {icon: 'clock', className: 'light'},
  [STATUS.PENDING]: {icon: 'clock', className: 'light'},
  [STATUS.ERROR]: {icon: 'times-circle', className: 'light red'}
};

const StepName = ({
  step,
  showToggle,
  open
}: {
  step: WorkflowStep | GroupedInputStep;
  showToggle?: boolean;
  open?: boolean;
}) => {
  let {state} = useContext(VulcanContext);
  const {workflow, status} = state;
  if (!workflow) return null;
  const statusStr = statusStringOfStepOrGroupedStep(step, workflow, status);
  let icon = icons[statusStr] || icons[STATUS.PENDING];

  let className = `step-status-icon ${icon.className}`;
  let IconComponent = (
    <Icon
      title={step.label || step.name}
      className={className}
      icon={icon.icon}
    />
  );

  if (STATUS.RUNNING === statusStr) {
    IconComponent = <AnimatedClock />;
  }

  return (
    <div className='step-name'>
      <div className='step-status-icon-wrapper'>{IconComponent}</div>
      <div className='step-button'>{step.label || step.name}</div>
      {showToggle ? open ? <ExpandLessIcon /> : <ExpandMoreIcon /> : null}
    </div>
  );
};

export default StepName;
