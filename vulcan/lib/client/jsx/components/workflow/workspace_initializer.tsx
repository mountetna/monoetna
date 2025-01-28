import React, {useCallback, useContext, useEffect, useState} from 'react';
import * as _ from 'lodash';

import {VulcanContext} from '../../contexts/vulcan_context';
import {
  uiContentsFromFiles,
  paramValuesFromRaw,
  vulcanConfigFromRaw,
  workflowByName,
  allFilesToBuffer,
  filesReturnToMultiFileContent,
  updateStepStatusesFromRunStatus,
  workspaceFromRaw
} from '../../selectors/workflow_selectors';

import WorkspaceManager from './session/workspace_manager';
import StepsList from './steps/steps_list';
import { paramValuesToRaw } from '../../selectors/workflow_selectors';
import { defaultWorkspaceStatus, FileContentResponse, Workspace, WorkspaceRaw, WorkspaceStatus } from '../../api_types';
import {
  setWorkflow,
  setAutoPassStep,
  setFullWorkspaceState
} from '../../actions/vulcan_actions';
import {
  defaultStepStatus,
  MultiFileContentResponse,
} from '../../api_types';
import {runPromise, useAsyncCallback} from 'etna-js/utils/cancellable_helpers';

export default function WorkspaceInitializer({
  workflowName,
  workspaceId,
  projectName
}: {
  workflowName: string;
  projectName: string;
  workspaceId: number;
}) {
  const {
    state,
    dispatch,
    showErrors,
    getWorkspace,
    getFileNames,
    readFiles
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

  const initializeFromWorkspace = useAsyncCallback(function* () {
    const workflow = workflowByName(workflowName, state);
    if (!!workflow) dispatch(setWorkflow(workflow, projectName));
    
    // workspace 
    const workspaceRaw: WorkspaceRaw = yield* runPromise(showErrors(getWorkspace(projectName, workspaceId)));
    if (!('workspace_id' in workspaceRaw)) {
      console.log("workspaceRaw is not a workspace");
    }
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

    // File pulls & inputUIs
    const fileNames: string[] = yield* runPromise(showErrors(getFileNames(projectName, workspaceId)));
    if (!Array.isArray(fileNames)) {
      console.log("fileNames is not an array");
    }
    status['output_files'] = fileNames;
    const fileGrabs = allFilesToBuffer(workspace).filter(f => fileNames.includes(f));
    if (fileGrabs.length > 0) {
      const filesContentRaw: MultiFileContentResponse = yield* runPromise(readFiles(projectName, workspaceId, fileGrabs))
      if (!Array.isArray(filesContentRaw)) {
        console.log("filesContentRaw is not an array");
      }
      const filesContent = filesReturnToMultiFileContent(filesContentRaw);
      status['file_contents'] = filesContent;
      status['ui_contents'] = uiContentsFromFiles(workspace, filesContent);
    } else {
      status['ui_contents'] = uiContentsFromFiles(workspace);
    }

    dispatch(setFullWorkspaceState(workspace, status));

    // Auto-pass for fully-defaulted params
    if (workspace.vignette?.includes('Primary inputs are skippable') && !workspace.last_job_status) {
      dispatch(setAutoPassStep(null));
    }
  }, [workflowName, state, projectName, dispatch]);

  useEffect(() => {
    if (!state.workspace) {
      // getLocalSession(workspaceId, projectName).then(
      // (localSession) => {
        // cancelPolling();

        // if (!!localSession) {
        //   initializeFromWorkspaceAndLocal(localSession);
        // } else {
          initializeFromWorkspace();
        // }
      // });
    }
  }, []);

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
