import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';
import {defaultSession} from '../../reducers/vulcan_reducer';
import {workflowByName} from '../../selectors/workflow_selectors';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';
import {
  setSession,
  setWorkflow,
  setSessionAndFigure
} from '../../actions/vulcan_actions';
import {json_get} from 'etna-js/utils/fetch';

export default function WorkflowManager({
  workflowName,
  figureId,
  projectName
}: {
  workflowName: string;
  projectName: string;
  figureId: number;
}) {
  const {
    state,
    dispatch,
    getLocalSession,
    cancelPolling,
    requestPoll
  } = useContext(VulcanContext);

  const {workflows} = state;

  useEffect(() => {
    if (workflows.length && projectName) {
      getLocalSession(workflowName, projectName, figureId).then((session) => {
        let workflow = workflowByName(workflowName, state);
        cancelPolling();

        if (!workflow) return;

        if (session) {
          dispatch(setWorkflow(workflow, projectName));
          dispatch(setSession(session));
        } else {
          if (figureId) {
            json_get(`/api/${projectName}/figure/${figureId}`).then(
              (figureResponse) => {
                if (workflow) dispatch(setWorkflow(workflow, projectName));
                dispatch(setSessionAndFigure(figureResponse));
              }
            );
          } else {
            dispatch(setWorkflow(workflow, projectName));
            let session = {
              ...defaultSession,
              workflow_name: `${workflowName}.cwl`,
              project_name: projectName
            };
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
      <SessionManager />
      <StepsList />
    </div>
  );
}
