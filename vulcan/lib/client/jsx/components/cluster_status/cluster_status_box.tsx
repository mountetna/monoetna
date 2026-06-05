import React, {useCallback, useState, useContext, useEffect} from 'react';

import Grid from '@material-ui/core/Grid';
import LatencyCheck from './latency_check';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  status_box: {
    width: 'auto',
    border: '1px solid #ccc',
    borderRadius: '2px',
    padding: '5px',
    backgroundColor: '#eee',
    boxShadow: 'inset 0px 0px 4px #ccc'
  }
}));

const ClusterStatusBox = ({children}:{
  children: any;
}) => {
  const classes = useStyles();

  return <Grid container alignItems='center'  className={ classes.status_box }>
    { children }
  </Grid>
}

export default ClusterStatusBox;
