import React, {useCallback, useState, useContext} from 'react';

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

  let {showErrors,
    showError,
    getConnectionLatency
  } = useContext(VulcanContext);

  const handleCheck = useCallback( () => {
    // Check if workflow name is available
    setChecking(true);
    showErrors(getConnectionLatency(projectName), () => {setChecking(false)})
    .then( (latencyReturn) => {
      setChecking(false);
      const latencyString = latencyReturn.latency;
      setLatency(parseFloat(latencyString.replace('ns', '')))
    })
  }, [projectName, getConnectionLatency]);

  console.log({latency})

  const valShow = 'Latency: ' + (latency!=null ? latency + 'ms' : '???ms');
  const actionIcon = checking ? <LoadingIcon/> : <TimerIcon/>;
  const tooltipText = (checking ? 'Measuring' : latency!=null ? 'Re-measure' : 'Measure') + ' latency in calculation server connection';

  return <Tooltip title={tooltipText}>
    <Button
      onClick={handleCheck}
      color='primary'
      variant='contained'
    >
      {valShow}
      {actionIcon}
    </Button>
  </Tooltip>
}

