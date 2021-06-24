import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';
import {defaultInputValues, workflowByName} from '../../selectors/workflow_selectors';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';
import {setInputs, setSession, setWorkflow} from "../../actions/vulcan_actions";
import {defaultVulcanSession} from "../../api_types";

export default function WorkflowManager({workflowName, projectName}: { workflowName: string, projectName: string }) {
  const {
    state,
    dispatch,
    getLocalSession,
    cancelPolling,
  } = useContext(VulcanContext);

  const workflow = workflowByName(workflowName, state);

  useEffect(() => {
    if (workflow && projectName) {
      dispatch(setWorkflow(workflow, projectName));

      getLocalSession(workflow, projectName).then((session) => {
        cancelPolling();

        if (!session) {
          // Set the default input values
          dispatch(setInputs(defaultInputValues(workflow)));
        } else {
          dispatch(setSession(session));
        }
      });
    }
  }, [dispatch, getLocalSession, projectName, workflow]);

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
