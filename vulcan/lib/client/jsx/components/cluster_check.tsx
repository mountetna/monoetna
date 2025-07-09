import React, {useCallback, useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../contexts/vulcan_context';
import Alert from '@material-ui/lab/Alert';
import AlertTitle from '@material-ui/lab/AlertTitle';

export default function ClusterStatusReport({projectName}: {
  projectName: string;
}) {
  const [connected, setConnected] = useState(true);
  const [message, setMessage] = useState('');
  const [expectedDown, setExpectedDown] = useState(false);

  let {showErrors,
    getClusterStatus
  } = useContext(VulcanContext);

  const handleCheck = useCallback( () => {
    showErrors(getClusterStatus(projectName))
    .then( (StatusReturn) => {
      setConnected(StatusReturn.connection_success)
      setExpectedDown(StatusReturn.expected_down)
      setMessage(StatusReturn.message)
    })
  }, [projectName, getClusterStatus]);

  // Run once, automatically on page load
  useEffect(() => {
    handleCheck()
  }, [])

  // Connected with no message to show
  if (message==='' && connected) return null

  // Connected, but message to show (likely: planned downage on the horizon)
  let level: 'info' | 'error' | 'warning' | 'success' = 'info';
  let highlight: string | undefined;
  let messageShow = <>{message}</>;
  // Disconnected, but expectedly
  if (!connected && expectedDown) {
    level = 'warning';
    highlight = 'Expected Connection Failure. Workspaces are inaccessible.';
  }
  // Disconnected and unexpectedly so
  if (!connected && !expectedDown) {
    level = 'error';
    highlight = 'Unexpected Connection Failure. Workspaces are inaccessible.';
    messageShow = <>Try again in a few minutes, but let the Data Library team know if this issue persists.<p>{'Also note: '+message}</p></>
  }

  return <Alert severity={level} style={{maxHeight: '44px', overflowY: 'auto'}}>
    {highlight && <AlertTitle>
      {highlight}
    </AlertTitle>}
    {messageShow}
  </Alert>
}