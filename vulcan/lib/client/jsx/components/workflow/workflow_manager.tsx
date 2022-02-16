import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../contexts/vulcan_context';
import {defaultSession} from '../../reducers/vulcan_reducer';
import {cwlName, workflowByName} from '../../selectors/workflow_selectors';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';
import {
  setSession,
  setWorkflow,
  setSessionAndFigure
} from '../../actions/vulcan_actions';
import {json_get} from 'etna-js/utils/fetch';
import {BufferedInputsContext} from '../../contexts/input_state_management';

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
  const {setInputs} = useContext(BufferedInputsContext);

  const {workflows} = state;

  useEffect(() => {
    if (workflows.length && projectName) {
      getLocalSession(workflowName, projectName, figureId).then((session) => {
        let workflow = workflowByName(workflowName, state);
        cancelPolling();

        if (session) {
          if (workflow) dispatch(setWorkflow(workflow, projectName));
          dispatch(setSession(session));
        } else if (figureId) {
          json_get(`/api/${projectName}/figure/${figureId}`).then(
            (figureResponse) => {
              workflow = workflowByName(figureResponse.workflow_name, state);
              if (workflow) dispatch(setWorkflow(workflow, projectName));
              dispatch(setSessionAndFigure(figureResponse));
              console.log('calling set inputs for buffer');
              // setInputs((inputs) => ({...inputs, ...figureResponse.inputs}));
            }
          );
        } else {
          if (workflow) dispatch(setWorkflow(workflow, projectName));
          let session = {
            ...defaultSession,
            workflow_name: cwlName(workflowName) || workflowName,
            project_name: projectName
          };
          dispatch(setSession(session));
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
