import React, {useCallback, useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../contexts/vulcan_context';
import Alert from '@material-ui/lab/Alert'

export default function ClusterKnownStatusReport({projectName}: {
  projectName: string;
}) {
  const [message, setMessage] = useState('');
  const [expectedDown, setExpectedDown] = useState(false);

  let {showErrors,
    getClusterStatus
  } = useContext(VulcanContext);

  const handleCheck = useCallback( () => {
    showErrors(getClusterStatus(projectName))
    .then( (StatusReturn) => {
      setExpectedDown(StatusReturn.expected_down)
      setMessage(StatusReturn.message)
    })
  }, [projectName, getClusterStatus]);

  useEffect(() => {
    handleCheck()
  }, [])

  if (message==='') return null
  return <Alert severity={expectedDown ? 'warning' : 'info'}>
    {message}
  </Alert>
}