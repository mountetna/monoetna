import React, {
  useContext, useState, useEffect, useCallback, useMemo
} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import StepIcon from './step_icon';
import {labelOfStepOrGroupedStep} from '../../../selectors/workflow_selectors';
import { statusOfStep } from '../../../selectors/workflow_selectors';
import StepUserInputDrawer from './step_user_input_drawer';
import {STATUS} from '../../../api_types';
import {WorkflowStepGroup} from "../user_interactions/inputs/input_types";

import Card from '@material-ui/core/Card';
import Typography from '@material-ui/core/Typography';
import Collapse from '@material-ui/core/Collapse';
import IconButton from '@material-ui/core/IconButton';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';

import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ExpandLessIcon from '@material-ui/icons/ExpandLess';

const useStyles = makeStyles((theme) => ({
  card: {
    borderRadius: 0,
    border: '1px solid #eee'
  },
  error: {
    border: '1px solid red'
  },
  label: {
    paddingLeft: '5px'
  },
  header: {
    padding: '5px 5px 5px 10px'
  }
}));

export default function StepUserInputWrapper({group}: { group: WorkflowStepGroup }) {
  const [open, setOpen] = useState(true);
  const {state} = useContext(VulcanContext);

  const {status} = state;

  const toggleInputs = useCallback(() => setOpen(!open), [setOpen, open]);
  const allInnerStatus = useMemo(() => group.steps.map((step) => statusOfStep(step, status)), [group, status]);
  const allStepsComplete = useMemo(() => allInnerStatus.every((s) => s && s.status === STATUS.COMPLETE),
    [allInnerStatus]
  );

  const hasValidationErrors = group.steps.some((step) => state.validationErrors.some(([stepName]) => stepName === step.name));
  const shouldOpen = hasValidationErrors;
  const shouldClose = allStepsComplete && !hasValidationErrors;

  useEffect(() => {
    if (shouldClose) {
      setOpen(false);
    } else if (shouldOpen) {
      setOpen(true);
    }
  }, [shouldOpen, shouldClose]);

  const classes = useStyles();

  const label = labelOfStepOrGroupedStep(group);
  return (<Card elevation={0} className={`${classes.card} ${hasValidationErrors ? classes.error : '' } step-user-input`}>
    <Grid justify='space-between' container onClick={toggleInputs}>
      <Grid item container style={{width:'auto'}}>
        <StepIcon step={group}/>
        <Typography className={classes.label}>{label}</Typography>
      </Grid>
      <IconButton size='small' onClick={ () => setOpen(!open) }>
        { open ? <ExpandLessIcon fontSize='small'/> : <ExpandMoreIcon fontSize='small'/>}
      </IconButton>
    </Grid>
    <Collapse className='step-user-input-inputs' in={open}>
      <StepUserInputDrawer group={group}/>
    </Collapse>
  </Card>);
}
