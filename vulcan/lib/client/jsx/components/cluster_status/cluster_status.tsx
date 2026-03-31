import React, {useCallback, useState, useContext, useEffect} from 'react';

import Grid from '@material-ui/core/Grid';
import LatencyCheck from './latency_check';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  status: {
    width: 'auto',
    marginLeft: '25px'
  }
}));

const ClusterStatus = ({projectName}: {
  projectName: string;
}) => {
  const classes = useStyles();

  return <Grid container alignItems='center'  className={ classes.status }>
    <LatencyCheck/>
  </Grid>
}

export default ClusterStatus;
