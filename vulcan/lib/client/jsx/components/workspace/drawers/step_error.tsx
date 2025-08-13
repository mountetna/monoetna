import React, {useContext, useState} from 'react';
import {VulcanContext} from '../../../contexts/vulcan_context';
import {statusOfStep} from '../../../selectors/workflow_selectors';
import {WorkspaceStep} from '../../../api_types';
import { useInputFeedStyles } from '../input_feed';
import Card from '@material-ui/core/Card';
import Grid from '@material-ui/core/Grid';
import StepIcon from './step_elements/step_icon';
import Typography from '@material-ui/core/Typography';
import IconButton from '@material-ui/core/IconButton';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ExpandLessIcon from '@material-ui/icons/ExpandLess';
import Collapse from '@material-ui/core/Collapse';

export default function StepError({step}: {step: WorkspaceStep}) {
  let {state} = useContext(VulcanContext);
  const stepStatus = statusOfStep(step, state.status, state.workspace);

  const [open, setOpen] = useState(true);

  const classes = useInputFeedStyles();

  if (!stepStatus) return null;

  let message = stepStatus.error || 'Something went wrong with this step.  We haven\'t automated log retrieval yet, so reach out to Dan for help.';

  return (
    <Card
      elevation={0}
      className={`${classes.card} ${classes.error} step-error`}
    >
      <Grid
        className={classes.header}
        justifyContent='space-between'
        container
        onClick={() => setOpen(!open)}
      >
        <Grid item container style={{width: 'auto'}}>
          <StepIcon step={step} />
          <Typography className={classes.label}>{step.name}</Typography>
        </Grid>
        <IconButton size='small' onClick={() => setOpen(!open)}>
          {open ? (
            <ExpandLessIcon fontSize='small' />
          ) : (
            <ExpandMoreIcon fontSize='small' />
          )}
        </IconButton>
      </Grid>
      <Collapse className='step-user-input-inputs' in={open}>
        <Grid style={{paddingLeft: '10px', paddingRight: '10px'}}>
          <Typography>{'Error:'}</Typography>
          <textarea readOnly value={message}/>
        </Grid>
      </Collapse>
    </Card>
  )
}
