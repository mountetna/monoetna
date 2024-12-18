import React, {useCallback, useContext, useEffect, useState} from 'react';
import * as _ from 'lodash';

import {VulcanContext} from '../../contexts/vulcan_context';
import {
  uiContentsFromFiles,
  paramValuesFromRaw,
  vulcanConfigFromRaw,
  workflowByName,
  allFilesToBuffer,
  filesReturnToMultiFileContent
} from '../../selectors/workflow_selectors';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';
import { paramValuesToRaw } from '../../selectors/workflow_selectors';
import { FileContentResponse } from '../../api_types';
import {
  setWorkflow,
  setAutoPassStep,
  setWorkspace,
  setStatusFromStatuses,
  setUIValues,
  setLastConfig,
  setWorkspaceFiles,
  setFileContent,
  setFilesContent
} from '../../actions/vulcan_actions';
import {
  defaultStepStatus,
  MultiFileContent,
  MultiFileContentResponse,
  WorkspaceRaw
} from '../../api_types';

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
    getLocalSession,
    cancelPolling,
    requestPoll,
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

  const initializeFromWorkspace = useCallback(() => {
    const workflow = workflowByName(workflowName, state);
    if (!!workflow) dispatch(setWorkflow(workflow, projectName));
    showErrors(getWorkspace(projectName, workspaceId).then(
      (workspaceRaw: WorkspaceRaw | Response) => {
        // ToDo: Fix this error handling!!!!!
        if (!('workspace_id' in workspaceRaw)) return;
        const vulcanConfig = vulcanConfigFromRaw(workspaceRaw.vulcan_config);
        const workspace = {
          ...workspaceRaw,
          vulcan_config: vulcanConfig
        };
        dispatch(setWorkspace(workspace, projectName));
        dispatch(setStatusFromStatuses(
          !!workspace.last_job_status ?
            workspace.last_job_status :
            Object.fromEntries(workspace.dag.map(
              stepName => [stepName, defaultStepStatus.statusFine]
            ))
        ));
        const param_vals = paramValuesFromRaw(workspace.last_config, workspace);
        dispatch(setUIValues(param_vals));
        dispatch(setLastConfig(
          !!workspace.last_config ?
            workspace.last_config :
            paramValuesToRaw(param_vals)
        ))
        getFileNames(projectName, workspaceId).then(
          (fileNames: string[] | Response) => {
            if (!Array.isArray(fileNames)) return;
            dispatch(setWorkspaceFiles(fileNames));
            const fileGrabs = allFilesToBuffer(workspace).filter(f => fileNames.includes(f));
            if (fileGrabs.length > 0) {
              readFiles(projectName, workspaceId, fileGrabs).then(
                (filesContentRaw: MultiFileContentResponse | Response) => {
                  if (!Array.isArray(filesContentRaw)) return;
                  const filesContent = filesReturnToMultiFileContent(filesContentRaw);
                  dispatch(setFilesContent(filesContent))
                  dispatch(setUIValues(uiContentsFromFiles(workspace, filesContent)));
              })
            } else {
              dispatch(setUIValues(uiContentsFromFiles(workspace)));
            }
        })
        if (workflow && workflow.vignette?.includes('Primary inputs are skippable') && !workspace.last_job_status) {
          dispatch(setAutoPassStep(null));
        }
      })
    )
  }, [workflowName, state, projectName, dispatch]);

  useEffect(() => {
    if (!state.workspace) {
      getLocalSession(workspaceId, projectName).then(
      (localSession) => {
        // cancelPolling();

        // if (!!localSession) {
        //   initializeFromWorkspaceAndLocal(localSession);
        // } else {
          initializeFromWorkspace();
        // }
      });
    }
  }, []);

  if (!state.workflow) {
    return null;
  }

  return (
    <div className='workspace-manager'>
      <SessionManager key={workspaceId} />
      <StepsList />
    </div>
  );
}
