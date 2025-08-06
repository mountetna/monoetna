import React, {useCallback, useState, useContext, useEffect} from 'react';

import Tooltip from '@material-ui/core/Tooltip';
import Button from '@material-ui/core/Button';
import {VulcanContext} from '../../contexts/vulcan_context';
import LoadingIcon from '../dashboard/loading_icon';
import TimerIcon from '@material-ui/icons/Timer';

export default function LatencyCheckButton({projectName}: {
  projectName: string;
}) {
  const [checking, setChecking] = useState(false);
  const [latency, setLatency] = useState<number | null>(null);

  let { getConnectionLatency } = useContext(VulcanContext);

  useEffect( () => {
    const getLatency = async () => {
      const { latency } = await getConnectionLatency(projectName);
      setLatency( Math.round(latency / 10) / 100 );
    }

    getLatency();
    const timer = setInterval( getLatency, 5000 );

    return () => clearInterval(timer);
  }, [] );

  const valShow = 'Latency: ' + (latency!=null ? latency + 's' : '???s');
  const actionIcon = checking ? <LoadingIcon/> : <TimerIcon/>;
  const tooltipText = (checking ? 'Measuring' : latency!=null ? 'Re-measure' : 'Measure') + ' latency in calculation server connection';

  return <Tooltip title={tooltipText}>
    <Button
      color='primary'
      variant='contained'
    >
      {valShow}
      {actionIcon}
    </Button>
  </Tooltip>
}
