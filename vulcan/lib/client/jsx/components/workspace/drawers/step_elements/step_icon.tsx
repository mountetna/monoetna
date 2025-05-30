import React, {useContext} from 'react';

import {VulcanContext} from '../../../../contexts/vulcan_context';
import AnimatedClock from './animated_clock';
import {WorkspaceStep} from '../../../../api_types';
import {statusStringOfStepOrGroupedStep} from '../../../../selectors/workflow_selectors';
import {STATUS} from '../../../../api_types';
import {WorkspaceStepGroup} from '../../ui_definitions/input_types';
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
    color: 'black'
  },
  [STATUS.PENDING]: {
    color: 'gray'
  },
  [STATUS.UPCOMING]: {
    color: 'black'
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
  [STATUS.UPCOMING]: {icon: ClockIcon},
  [STATUS.ERROR]: {icon: ErrorIcon}
};
const StepIcon = ({
  step
}: {
  step: WorkspaceStep | WorkspaceStepGroup;
}) => {
  let {state} = useContext(VulcanContext);
  const {workflow, status} = state;
  const classes = useStyles();
  if (!workflow) return null;
  const statusStr = statusStringOfStepOrGroupedStep(step, state.workspace, status);
  let icon_config = icons[statusStr] || icons[STATUS.PENDING];

  let IconComponent:IconClass = STATUS.RUNNING === statusStr ? AnimatedClock : icon_config.icon;

  return <Tooltip title={statusStr}>
    <IconComponent className={`${classes[ statusStr ]} ${classes.icon} icon`}/>
  </Tooltip>;
};

export default StepIcon;
