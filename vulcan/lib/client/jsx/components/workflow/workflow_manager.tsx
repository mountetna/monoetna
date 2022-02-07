import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';
import {workflowByName} from '../../selectors/workflow_selectors';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';
import {setSession, setWorkflow} from "../../actions/vulcan_actions";
import { json_get } from 'etna-js/utils/fetch';

export default function WorkflowManager({workflowName, figureId, projectName}: { workflowName: string, projectName: string }) {
  const {
    state,
    dispatch,
    getLocalSession,
    cancelPolling,
    requestPoll,
    setSession,
  } = useContext(VulcanContext);

  const [ currentWorkflowName, setCurrentWorkflowName ] = useState(workflowName)
  
  const workflow = currentWorkflowName ? workflowByName(currentWorkflowName.replace('.cwl',''), state) : undefined;

  useEffect(() => {
    if (figureId && !currentWorkflowName) {
      // this might fire too often
      json_get(`/api/${projectName}/figure/${figureId}`).then( figure => {
        setCurrentWorkflowName( figure.workflow_name ); setSession(figure);
      })
    }
  }, []);

  useEffect(() => {
    if (workflow && projectName) {
      getLocalSession(workflow, projectName).then((session) => {
        cancelPolling();

        dispatch(setWorkflow(workflow, projectName));

        if (session) {
          dispatch(setSession(session));
        }

        requestPoll();
      });
    }
  }, [cancelPolling, dispatch, getLocalSession, projectName, workflow]);

  if (!state.workflow) {
    return null;
  }

  return (
    <div className='workflow-manager'>
      <SessionManager/>
      <StepsList/>
    </div>
  );
}
