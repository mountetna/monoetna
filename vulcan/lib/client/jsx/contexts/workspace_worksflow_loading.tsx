import {Dispatch} from 'react';
import {setWorkflowsWorkspaces, VulcanAction} from '../actions/vulcan_actions';
import {defaultApiHelpers} from './api';
import { runPromise, useAsync } from 'etna-js/utils/cancellable_helpers';

export function useWorkspaceWorkflowLoading(
    key: string,
    dispatch: Dispatch<VulcanAction>,
    getWorkflows: typeof defaultApiHelpers.getWorkflows,
    getWorkspaces: typeof defaultApiHelpers.getWorkspaces,
    getWorkspace: typeof defaultApiHelpers.getWorkspace,
    showErrors: typeof defaultApiHelpers.showErrors,
    projectName: string | undefined,
    workspaceId: number | null,
    run_update: boolean
) {

    useAsync(function* () {
        if (run_update) {
            const update = {};
            if (!!projectName) {
                update['workspaces'] = (yield* runPromise(showErrors(getWorkspaces(projectName)))).workspaces
                update['workflows'] = (yield* runPromise(showErrors(getWorkflows(projectName)))).workflows
            }
            // We WANT to rely on the workspace_initializer for this part!
            // if (!!projectName && !!workspaceId) {
            //     update['workspace'] = (yield* runPromise(showErrors(getWorkspace(projectName, workspaceId))))
            // }
            if (Object.keys(update).length > 0) {
                dispatch(setWorkflowsWorkspaces(update));
            }
        }
    }, [key, dispatch, showErrors, getWorkspace, getWorkflows, getWorkspaces, projectName, workspaceId, run_update]);
}