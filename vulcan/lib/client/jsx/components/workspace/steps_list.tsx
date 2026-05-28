import React, {useContext, useState, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import {WorkspaceContext} from '../../contexts/workspace_context';
import {outputUINames, outputUIsWithInputsReady, stepOfName} from '../../selectors/workflow_selectors';
import {useWorkspace} from '../../contexts/workspace_context';
import { VulcanConfig } from '../../api_types';
import StepIconName from './drawers/step_elements/step_icon_name';
import FlatButton from 'etna-js/components/flat-button'
import Grid from '@material-ui/core/Grid';
import Drawer from '@material-ui/core/Drawer';
import IconButton from '@material-ui/core/IconButton';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import Typography from '@material-ui/core/Typography';
import ChevronRightIcon from '@material-ui/icons/ChevronRight';
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft';

const useStyles = makeStyles((theme) => ({
  steps: {
    borderLeft: '1px solid gray'
  },
  stepsOpen: {
    width: '300px'
  },
  stepsClosed: {
    width: '60px'
  }
}));

const Step = ({name}) => {
  return <Grid sx={{ height: '60px' }}>
  {name}
  </Grid>
}

export default function StepsList() {
  const [open, setOpen] = useState(true);
  const {state: { file_contents, workspaceDetails, dag } } = useContext(WorkspaceContext);

  const classes = useStyles();

  function handleToggle() {
    setOpen(!open);
  }

  console.log({workspaceDetails, dag});

  const hasCompletedOutputs = useMemo(() => {
    if (!workspaceDetails) return false;
    return outputUIsWithInputsReady(workspaceDetails, file_contents).length > 0
  }, [workspaceDetails, file_contents]);

  //useEffect(() => {
    //setOpen(!hasCompletedOutputs || state.isSyncing);
  //}, [state.isSyncing, hasCompletedOutputs]);

  if (!workspaceDetails || Object.keys(workspaceDetails).length == 0) return null;

  const stepNamesToStepIconNames = (step, index) => (
    <div key={index} className='step'>
      <StepIconName step={stepOfName(step, workspaceDetails.vulcan_config)}/>
    </div>
  )

  return (
    <Grid className={ `${classes.steps} ${open ? classes.stepsOpen : classes.stepsClosed}` } >
      <Grid 
        style={{ borderBottom: '1px solid gray' }}
        alignItems='center' container>
        <IconButton onClick={ () => setOpen(!open) }>
          { open ? <ChevronRightIcon/> : <ChevronLeftIcon/> }
        </IconButton>
        { open && <Typography component="span">Progress</Typography> }
      </Grid>
      <List>
        {workspaceDetails.dag_flattened.map(
          stepName => <ListItem key={stepName}>
            <Step name={stepName}/>
          </ListItem>
        )}
      </List>
    </Grid>
  );
};
