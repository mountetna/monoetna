import React, {useCallback, useState, useContext, useEffect} from 'react';

import Tooltip from '@material-ui/core/Tooltip';
import Grid from '@material-ui/core/Grid';
import {VulcanContext} from '../../contexts/vulcan_context';
import LoadingIcon from '../dashboard/loading_icon';
import TimerIcon from '@material-ui/icons/Timer';
import {makeStyles} from '@material-ui/core/styles';
import ClusterStatusBox from './cluster_status_box';
import {getConnectionLatency} from '../../contexts/api';

const useStyles = makeStyles((theme) => ({
  latency: {
    fontSize: '0.9em',
    color: theme.palette.secondary.main
  }
}));
export default function LatencyCheck({projectName}: {
  projectName: string;
}) {
  const [latency, setLatency] = useState<number | null>(null);

  const classes = useStyles();

  const getLatency = async () => {
    const { latency } = await getConnectionLatency();
    setLatency( Math.round(latency / 100) / 10 );
  }

  useEffect( () => {
    getLatency();
    const timer = setInterval( getLatency, 120000 );

    return () => clearInterval(timer);
  }, [] );

  const valShow = (latency != null ? latency + 's' : '???s');

  return <ClusterStatusBox>
    <Tooltip title='Cluster latency (delay)'>
      <Grid container alignItems='center' className={ classes.latency } onClick={() => getLatency()}>
        <TimerIcon fontSize="small"/>
        {valShow}
      </Grid>
    </Tooltip>
  </ClusterStatusBox>
}
