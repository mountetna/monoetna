import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import Icon from 'etna-js/components/icon';
import AnimatedClock from './animated_clock';
import {WorkflowStep} from '../../../api_types';
import {labelOfStepOrGroupedStep, statusStringOfStepOrGroupedStep} from '../../../selectors/workflow_selectors';
import {STATUS} from '../../../api_types';
import {WorkflowStepGroup} from "../user_interactions/inputs/input_types";
import {makeStyles} from '@material-ui/core/styles';
import Tooltip from '@material-ui/core/Tooltip';
import CheckIcon from '@material-ui/icons/Check';
import ClockIcon from '@material-ui/icons/AccessTime';
import ErrorIcon from '@material-ui/icons/HighlightOff';

const useStyles = makeStyles((theme) => ({
  icon: {
    opacity: 0.5
  },
  [STATUS.COMPLETE]: {
    color: 'green'
  },
  [STATUS.RUNNING]: {
    color: 'gray'
  },
  [STATUS.PENDING]: {
    color: 'gray'
  },
  [STATUS.ERROR]: {
    color: 'red'
  }
}));

type IconClass = typeof CheckIcon;

const icons: {[k: string]: {icon: IconClass}} = {
  [STATUS.COMPLETE]: {icon: CheckIcon},
  [STATUS.RUNNING]: {icon: ClockIcon},
  [STATUS.PENDING]: {icon: ClockIcon},
  [STATUS.ERROR]: {icon: ErrorIcon}
};
const StepIcon = ({
  step
}: {
  step: WorkflowStep | WorkflowStepGroup;
}) => {
  let {state} = useContext(VulcanContext);
  const {workflow, status} = state;
  if (!workflow) return null;
  const statusStr = statusStringOfStepOrGroupedStep(step, workflow, status);
  const label = labelOfStepOrGroupedStep(step);
  let icon_config = icons[statusStr] || icons[STATUS.PENDING];

  let IconComponent:IconClass = STATUS.RUNNING === statusStr ? AnimatedClock : icon_config.icon;

  const classes = useStyles();
  return <Tooltip title={label}>
    <IconComponent className={`${classes[ statusStr ]} ${classes.icon}`}/>
  </Tooltip>
}

export default StepIcon;
