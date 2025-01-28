import React, { useState, useCallback } from 'react';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';
import WorkflowPane, { WorkflowPaneHeader } from './workflow-pane';

const useStyles = makeStyles( theme => ({
  confirm: {
    justifyContent: 'space-between',
    width: 'auto',
  }
}));
const RemovePane = ({update,selected}:{update:Function,selected:string|null}) => {
  const classes = useStyles();

  return <WorkflowPane mode='remove' selected={selected}>
    <WorkflowPaneHeader title='Disable this loader?'>
      <Grid item container className={ classes.confirm }>
      <Button onClick={ () => update({ disabled: true }) }>Remove</Button>
      </Grid>
    </WorkflowPaneHeader>
  </WorkflowPane>;
};

export default RemovePane;
