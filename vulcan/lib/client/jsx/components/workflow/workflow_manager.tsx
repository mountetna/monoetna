import React, {useCallback, useContext, useEffect, useState} from 'react';
import * as _ from 'lodash';

import {VulcanContext} from '../../contexts/vulcan_context';
import {defaultSession} from '../../reducers/vulcan_reducer';
import {
  cwlName,
  defaultInputs,
  selectFigure,
  selectSession,
  workflowByName
} from '../../selectors/workflow_selectors';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';
import {
  setSession,
  setWorkflow,
  setSessionAndFigureSeparately
} from '../../actions/vulcan_actions';
import {
  defaultFigure,
  VulcanFigure,
  VulcanFigureSession,
  VulcanSession
} from '../../api_types';

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
    requestPoll,
    showErrors,
    fetchFigure
  } = useContext(VulcanContext);

  const initializeFromSessionAndFigure = useCallback(
    (session: VulcanSession, figure: VulcanFigure, post: boolean = false) => {
      const workflow = workflowByName(session.workflow_name, state);
      if (workflow) dispatch(setWorkflow(workflow, projectName));
      dispatch(setSessionAndFigureSeparately(figure, session));

      requestPoll(post);
    },
    [projectName, dispatch, state, requestPoll]
  );

  const initializeFromFigure = useCallback(
    (figureId: number, localSession: VulcanFigureSession | null) => {
      showErrors(
        fetchFigure(projectName, figureId).then((figureResponse) => {
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
            selectSession(
              !fromDatabase && localSession ? localSession : figureResponse
            ),
            selectFigure(figureResponse),
            true
          );
        })
      );
    },
    [projectName, showErrors, fetchFigure, initializeFromSessionAndFigure]
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

  const initializeFromStoredSession = useCallback(
    (localSession: VulcanFigureSession) => {
      const workflow = workflowByName(localSession.workflow_name, state);
      if (!workflow) {
        initializeNewSession();
        return;
      }

      dispatch(setWorkflow(workflow, projectName));
      const defaults = defaultInputs(workflow);

      let discardStoredSession = true;

      if (
        !_.isEqual(defaults, localSession.inputs) &&
        Object.keys(localSession.inputs) >= Object.keys(defaults)
      ) {
        discardStoredSession = confirm(
          'You have an edited, unsaved version of this workflow. Discard it?'
        );
      }

      if (discardStoredSession) {
        initializeNewSession();
      } else {
        initializeFromSessionAndFigure(
          selectSession(localSession),
          defaultFigure
        );
      }
    },
    [
      projectName,
      state,
      dispatch,
      initializeFromSessionAndFigure,
      initializeNewSession
    ]
  );

  useEffect(() => {
    getLocalSession(workflowName, projectName, figureId).then(
      (localSession) => {
        cancelPolling();

        if (localSession && !figureId) {
          initializeFromStoredSession(localSession);
        } else if (figureId) {
          initializeFromFigure(figureId, localSession);
        } else {
          initializeNewSession();
        }
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
