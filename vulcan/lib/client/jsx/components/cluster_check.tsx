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

  useEffect(() => {
    handleCheck()
  }, [])

  if (message==='' && connected) return null

  let connectionMessage = !connected ? (expectedDown ? 'Expected' : 'Unexpected') + ' Connection Failure. Workspaces are inaccessible.' :
    undefined;

  return <Alert severity={!connected && expectedDown ? 'warning' : connected ? 'info' : 'error'} style={{maxHeight: '44px', overflowY: 'auto'}}>
    {connectionMessage && <AlertTitle>
      {connectionMessage}
    </AlertTitle>}
    {message}
  </Alert>
}