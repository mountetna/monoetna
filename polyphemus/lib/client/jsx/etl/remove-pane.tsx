import React, { useState, useCallback } from 'react';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';
import EtlPane, { EtlPaneHeader } from './etl-pane';

const useStyles = makeStyles( theme => ({
  confirm: {
    justifyContent: 'space-between',
    width: 'auto',
  }
}));
const RemovePane = ({update,selected}:{update:Function,selected:string|null}) => {
  const classes = useStyles();

  return <EtlPane mode='remove' selected={selected}>
    <EtlPaneHeader title='Disable this loader?'>
      <Grid item container className={ classes.confirm }>
        <Button onClick={ () => update({ archived: true }) }>Remove</Button>
      </Grid>
    </EtlPaneHeader>
  </EtlPane>
}

export default RemovePane;
