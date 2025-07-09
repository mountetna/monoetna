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
    showError,
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

  console.log({
    connected,
    expectedDown,
    message
  })

  // Connected with no message to show
  if (message==='' && connected) return null

  // Connected, but message to show (likely: planned downage on the horizon)
  let messageUse: string = message;
  // Disconnected, but expectedly
  if (!connected && expectedDown) {
    messageUse = `Expected Vulcan Connection Failure. Workspaces are inaccessible.\n${message}`;
  }
  // Disconnected and unexpectedly so
  if (!connected && !expectedDown) {
    messageUse = 'Unexpected Vulcan Connection Failure. Workspaces are inaccessible.\n\n' +
      `Try again in a few minutes, but let the Data Library team know if the issue persists.\n\nAlso note: ${message}`
  }
  showError(messageUse);
  return null
}