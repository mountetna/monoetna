import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';
import {workflowByName} from '../../selectors/workflow_selectors';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';
import {setSession, setWorkflow} from "../../actions/vulcan_actions";

export default function WorkflowManager({workflowName, projectName}: { workflowName: string, projectName: string }) {
  const {
    state,
    dispatch,
    getLocalSession,
    cancelPolling,
    requestPoll,
  } = useContext(VulcanContext);

  const workflow = workflowByName(workflowName, state);

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
