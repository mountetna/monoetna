import React, {useCallback, useContext, useEffect, useState} from 'react';
import * as _ from 'lodash';

import {VulcanContext} from '../../contexts/vulcan_context';
import {
  paramValuesFromRaw,
  workflowByIdFromWorkflows,
  paramValuesToRaw,
  updateStepStatusesFromRunStatus,
  workspaceFromRaw
} from '../../selectors/workflow_selectors';

import WorkspaceManager from './workspace_manager';
import StepsList from './steps_list';
import { defaultWorkspaceStatus } from '../../api_types';
import {
  setWorkflow,
  setAutoPassStep,
  setFullWorkspaceState
} from '../../actions/vulcan_actions';
import {
  defaultStepStatus,
} from '../../api_types';
import {runPromise, useAsyncCallback} from 'etna-js/utils/cancellable_helpers';

export default function WorkspaceInitializer({
  workspaceId,
  projectName
}: {
  projectName: string;
  workspaceId: number;
}) {
  const {
    state,
    dispatch,
    showErrors,
    getWorkspace,
  } = useContext(VulcanContext);

  // For now, we will ALWAYS initialize form the state of the workspace on c4
  // ToDo: Decide workings of using local version instead & implement.

  const useLocalPrompt =
    "You have unsaved changes. Click 'OK' to use your local version, or 'Cancel' to discard all unsaved changes and load the last saved state.";

  // const initializeFromWorkspaceAndLocal = useCallback(
  //   (workspaceId: number, localSession: VulcanFigureSession | null) => {
  //     showErrors(
  //       fetchFigure(projectName, figureId).then((figureResponse) => {
  //         let useLocal = true;

  //         if (
  //           localSession &&
  //           !_.isEqual(localSession.inputs, figureResponse.inputs)
  //         ) {
  //           useLocal = confirm(useLocalPrompt);
  //         }

  //         initializeFromSessionAndFigure(
  //           selectSession(
  //             useLocal && localSession ? localSession : figureResponse
  //           ),
  //           selectFigure(figureResponse)
  //         );
  //       })
  //     );
  //   },
  //   [projectName, showErrors, fetchFigure, initializeFromSessionAndFigure]
  // );

  const [initializeFromWorkspace] = useAsyncCallback(function* () {
    // workspace
    showErrors(getWorkspace(projectName, workspaceId))
    .then((workspaceRaw) => {
      const workspace = workspaceFromRaw(workspaceRaw);

      const status = defaultWorkspaceStatus;
      
      // step statuses
      const defaultStepStatuses = Object.fromEntries(workspace.dag.map(
        stepName => [stepName, defaultStepStatus]
      ))
      status['steps'] = !!workspace.last_job_status ?
        updateStepStatusesFromRunStatus(workspace.last_job_status, defaultStepStatuses).newStepStatus :
        defaultStepStatuses;

      // paramUIs
      const param_vals = paramValuesFromRaw(workspace.last_config, workspace);
      status['last_params'] = !!workspace.last_config ?
        workspace.last_config :
        paramValuesToRaw(param_vals)
      status['params'] = param_vals
      
      // Send it, with the true here triggering files to be updated in a next render
      dispatch(setFullWorkspaceState(workspace, status, true));

      // // Auto-pass for fully-defaulted params
      // if (workspace.vignette?.includes('Primary inputs are skippable') && !workspace.last_job_status) {
      //   dispatch(setAutoPassStep(Object.keys(param_vals)));
      // }
    })
  }, [projectName, workspaceId, dispatch]);

  useEffect(() => {
    if (state.workspace==null) {
      // getLocalSession(workspaceId, projectName).then(
      // (localSession) => {
        // cancelPolling();

        // if (!!localSession) {
        //   initializeFromWorkspaceAndLocal(localSession);
        // } else {
          initializeFromWorkspace();
        // }
      // });
    } else if (state.workflow.name == '') {
      const workflow = workflowByIdFromWorkflows(state.workspace.workflow_id, state.workflows);
      if (!!workflow) dispatch(setWorkflow(workflow, projectName));
    }
  }, [state.workspace, state.workflow.name]);

  if (!state.workflow) {
    return null;
  }

  return (
    <div className='workspace-manager'>
      <WorkspaceManager key={workspaceId} />
      <StepsList />
    </div>
  );
}
