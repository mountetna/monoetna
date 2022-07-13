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

  const selectWorkflowSnapshot = useCallback(
    (
      figure: VulcanFigure | VulcanFigureSession,
      defaultWorkflowName: string
    ) => {
      return figure.workflow_snapshot
        ? figure.workflow_snapshot
        : workflowByName(defaultWorkflowName, state);
    },
    [state]
  );

  const initializeFromSessionAndFigure = useCallback(
    (session: VulcanSession, figure: VulcanFigure) => {
      const workflow = selectWorkflowSnapshot(figure, session.workflow_name);
      if (workflow) dispatch(setWorkflow(workflow, projectName));
      dispatch(setSessionAndFigureSeparately(figure, session));

      requestPoll();
    },
    [projectName, dispatch, selectWorkflowSnapshot, requestPoll]
  );

  const useLocalPrompt =
    "You have unsaved changes. Click 'OK' to use your local version, or 'Cancel' to discard all unsaved changes and load the last saved state.";

  const initializeFromFigure = useCallback(
    (figureId: number, localSession: VulcanFigureSession | null) => {
      showErrors(
        fetchFigure(projectName, figureId).then((figureResponse) => {
          let useLocal = true;

          if (
            localSession &&
            !_.isEqual(localSession.inputs, figureResponse.inputs)
          ) {
            useLocal = confirm(useLocalPrompt);
          }

          initializeFromSessionAndFigure(
            selectSession(
              useLocal && localSession ? localSession : figureResponse
            ),
            selectFigure(figureResponse)
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
      project_name: projectName,
      inputs: workflow ? defaultInputs(workflow) : defaultSession.inputs
    };
    dispatch(setSession(session));
  }, [workflowName, state, projectName, dispatch]);

  const initializeFromStoredSession = useCallback(
    (localSession: VulcanFigureSession) => {
      const workflow = selectWorkflowSnapshot(
        localSession,
        localSession.workflow_name
      );
      if (!workflow) {
        initializeNewSession();
        return;
      }

      dispatch(setWorkflow(workflow, projectName));
      const defaults = defaultInputs(workflow);

      let useLocal = true;

      if (
        !_.isEqual(defaults, localSession.inputs) &&
        Object.keys(localSession.inputs) >= Object.keys(defaults)
      ) {
        useLocal = confirm(useLocalPrompt);
      }

      if (!useLocal) {
        initializeNewSession();
      } else {
        initializeFromSessionAndFigure(
          selectSession(localSession),
          selectFigure(localSession)
        );
      }
    },
    [
      projectName,
      selectWorkflowSnapshot,
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
