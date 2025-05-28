import {Dispatch} from 'react';
import {setWorkflowsWorkspaces, VulcanAction} from '../actions/vulcan_actions';
import {defaultApiHelpers} from './api';
import { runPromise, useAsync } from 'etna-js/utils/cancellable_helpers';

export function useWorkspacesWorkflowLoading(
    run_update: boolean,
    dispatch: Dispatch<VulcanAction>,
    getWorkflows: typeof defaultApiHelpers.getWorkflows,
    getWorkspaces: typeof defaultApiHelpers.getWorkspaces,
    showErrors: typeof defaultApiHelpers.showErrors,
    projectName: string | undefined,
) {

    useAsync(function* () {
        if (run_update) {
            const update = {};
            if (!!projectName) {
                showErrors(getWorkspaces(projectName))
                .then((workspacesReturn) => {
                    update['workspaces'] = workspacesReturn.workspaces;
                    showErrors(getWorkflows(projectName))
                    .then((workflowsReturn) => {
                        update['workflows'] = workflowsReturn.workflows
                        dispatch(setWorkflowsWorkspaces(update));
                    })
                })
            }
        }
    }, [run_update, dispatch, getWorkflows, getWorkspaces, showErrors, projectName]);
}