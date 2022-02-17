import React, {useCallback, useContext, useEffect, useState} from 'react';
import * as _ from 'lodash';

import {VulcanContext} from '../../contexts/vulcan_context';
import {defaultSession} from '../../reducers/vulcan_reducer';
import {
  cwlName,
  selectFigure,
  selectSession,
  workflowByName
} from '../../selectors/workflow_selectors';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';
import {
  setSession,
  setWorkflow,
  setSessionAndFigure,
  setSessionAndFigureSeparately
} from '../../actions/vulcan_actions';
import {json_get} from 'etna-js/utils/fetch';
import {
  VulcanFigure,
  VulcanFigureSession,
  VulcanSession
} from '../../api_types';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

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
  const invoke = useActionInvoker();

  const handleErrorResponse = useCallback(
    (err: Promise<any>) => {
      err.then((e: any) => invoke(showMessages(e.errors || [e.error])));
    },
    [invoke]
  );

  const initializeFromSessionAndFigure = useCallback(
    (session: VulcanSession, figure: VulcanFigure) => {
      const workflow = workflowByName(session.workflow_name, state);
      if (workflow) dispatch(setWorkflow(workflow, projectName));
      dispatch(setSessionAndFigureSeparately(figure, session));

      requestPoll();
    },
    [projectName, dispatch, state, requestPoll]
  );

  const initializeFromFigure = useCallback(
    (figureId: number, localSession: VulcanFigureSession | null) => {
      json_get(`/api/${projectName}/figure/${figureId}`)
        .then((figureResponse) => {
          let fromDatabase = true;

          if (
            localSession &&
            !_.isEqual(localSession.inputs, figureResponse.inputs)
          ) {
            fromDatabase = confirm(
              'You have an edited, local version of that figure. Discard it?'
            );
          }
          initializeFromSessionAndFigure(
            selectSession(fromDatabase ? figureResponse : localSession),
            selectFigure(figureResponse)
          );
        })
        .catch(handleErrorResponse);
    },
    [projectName, handleErrorResponse, initializeFromSessionAndFigure]
  );

  const initializeNewSession = useCallback(() => {
    const workflow = workflowByName(workflowName, state);
    if (workflow) dispatch(setWorkflow(workflow, projectName));
    let session = {
      ...defaultSession,
      workflow_name: cwlName(workflowName) || workflowName,
      project_name: projectName
    };
    dispatch(setSession(session));
  }, [workflowName, state, projectName, dispatch]);

  useEffect(() => {
    getLocalSession(workflowName, projectName, figureId).then(
      (localSession) => {
        cancelPolling();

        if (localSession && !figureId) {
          initializeFromSessionAndFigure(
            selectSession(localSession),
            selectFigure(localSession)
          );
        } else if (figureId) {
          initializeFromFigure(figureId, localSession);
        } else {
          initializeNewSession();
        }

        requestPoll();
      }
    );
  }, []);

  if (!state.workflow) {
    return null;
  }

  return (
    <div className='workflow-manager'>
      <SessionManager key={state.session.key} />
      <StepsList />
    </div>
  );
}
