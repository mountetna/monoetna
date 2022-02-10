import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';
import {defaultSession} from '../../reducers/vulcan_reducer';
import {localStorageKey} from '../../contexts/session_storage';
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
    requestPoll
  } = useContext(VulcanContext);

  const { workflows } = state;

  useEffect(() => {
    if (workflows.length && projectName) {
      getLocalSession(workflowName, figureId, projectName).then((session) => {
        let workflow;

        cancelPolling();

        if (session) {
          workflow = workflowByName(session.workflow_name, state);
          dispatch(setWorkflow(workflow, projectName));
          dispatch(setSession(session));
        } else {
          if (figureId) {
            json_get(`/api/${projectName}/figure/${figureId}`).then(
              figure => {
                figure.key = `${figure.project_name}/${figure.figure_id}`;
                workflow = workflowByName(figure.workflow_name, state);
                dispatch(setWorkflow(workflow, projectName));
                dispatch(setSession(figure));
              }
            )
          } else {
            workflow = workflowByName(`${workflowName}.cwl`, state);
            dispatch(setWorkflow(workflow, projectName));
            let session = {
              ...defaultSession,
              workflow_name: `${workflowName}.cwl`,
              project_name: projectName
            }
            console.log({workflow, session});
            dispatch(setSession(session));
          }
        }

        requestPoll();
      });
    }
  }, [cancelPolling, dispatch, getLocalSession, projectName, workflows]);

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
