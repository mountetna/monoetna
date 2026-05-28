import React, {useCallback, useContext, useEffect, useState} from 'react';
import * as _ from 'lodash';

import {VulcanContext} from '../../contexts/vulcan_context';
import {
  paramValuesFromRaw,
  workflowByIdFromWorkflows,
  paramValuesToRaw,
  updateStepStatusesFromRunStatus,
  workspaceFromRaw,
  updateStepStatusesFromJobsState
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
import { LoadingIconWithText } from '../dashboard/loading_icon';

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
    getState,
    getIsRunning
  } = useContext(VulcanContext);

  const [initializeFromWorkspace] = useAsyncCallback(function* () {
    // workspace
    showErrors(getWorkspace(projectName, workspaceId))
    .then((workspaceRaw) => {
      const workspace = workspaceFromRaw(workspaceRaw);

      console.log({workspace})

      const status = defaultWorkspaceStatus;

      // paramUIs
      const param_vals = paramValuesFromRaw(workspace.last_config, workspace);
      status['last_params'] = !!workspace.last_config ?
        workspace.last_config :
        paramValuesToRaw(param_vals)
      status['params'] = param_vals

      // Check if running
      showErrors(getIsRunning(projectName, workspaceId))
      .then((isRunningReturn) => {
        const isRunning = isRunningReturn['running']

        // step statuses
        // ToDo: Bring back use of RunStatus / error knowledge.
        const defaultStepStatuses = Object.fromEntries(workspace.dag_flattened.map(
          stepName => [stepName, defaultStepStatus]
        ))
        status['steps'] = defaultStepStatuses;
        if (!!workspace.last_config_id) {
          showErrors(getState(projectName, workspace.last_config_id))
          .then((stateReturn) => {
            status.last_file_accounting = stateReturn.files;
            status.last_jobs_accounting = stateReturn.jobs;
            status['steps'] = updateStepStatusesFromJobsState(stateReturn.jobs, status['steps'])

            // Send it, with the true here triggering files to be updated in a next render
            dispatch(setFullWorkspaceState(workspace, status, true, isRunning));
          })
        } else {
          dispatch(setFullWorkspaceState(workspace, status, true, isRunning));
        }
      })
    })
  }, [projectName, workspaceId, dispatch]);

  useEffect(() => {
    if (state.workspace==null) {
      initializeFromWorkspace();
    } else if (state.workflow.name == '') {
      const workflow = workflowByIdFromWorkflows(state.workspace.workflow_id, state.workflows);
      if (!!workflow) dispatch(setWorkflow(workflow, projectName));
    }
  }, [state.workspace, state.workflow.name]);

  if (!state.workspace) {
    return (
      <LoadingIconWithText text='Retrieving Workspace Context'/>
    );
  }

  return (
    <div className='workspace-manager'>
    </div>
  );
}
