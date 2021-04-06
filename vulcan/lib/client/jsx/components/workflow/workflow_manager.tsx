import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';
import {defaultInputValues, workflowByName} from '../../selectors/workflow_selectors';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';
import {setInputs, setSession, setWorkflow} from "../../actions/vulcan";
import {defaultVulcanSession} from "../../api_types";

export default function WorkflowManager({workflowName}: { workflowName: string }) {
  const {
    state,
    dispatch,
    getLocalSession
  } = useContext(VulcanContext);

  const workflow = workflowByName(workflowName, state);

  useEffect(() => {
    if (workflow) {
      dispatch(setWorkflow(workflow));

      getLocalSession(workflow).then((session) => {
        if (!session) {
          // Set the default input values
          dispatch(setInputs(defaultInputValues(workflow)));
        } else {
          dispatch(setSession(session));
        }
      });
    }
  }, [workflow]);

  return (
      <div className='workflow-manager'>
        <SessionManager/>
        <StepsList/>
      </div>
  );
}
