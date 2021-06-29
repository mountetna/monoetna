import {SessionStatusResponse, Workflow, WorkflowsResponse} from "../api_types";
import {InputValidator} from '../components/workflow/user_interactions/inputs/input_types';

function actionObject<T extends string, P>(type: T, payload: P): { type: T } & P {
    return { ...payload, type };
}

export function setWorkflows(workflows: WorkflowsResponse['workflows']) {
    return actionObject('SET_WORKFLOWS', {workflows});
}

export function setWorkflow(workflow: Workflow, projectName: string) {
    return actionObject('SET_WORKFLOW', {workflow, projectName});
}

export function setStatus(status: SessionStatusResponse['status']) {
    return actionObject('SET_STATUS', {status});
}

export function setDownloadedData(url: string, data: any) {
    return actionObject('SET_DOWNLOAD', {url, data});
}

export function releaseDownloadedData(url: string) {
    return actionObject('RELEASE_DOWNLOAD', {url});
}

export function setSession(session: SessionStatusResponse['session']) {
    return actionObject('SET_SESSION', {session});
}

// Fully replaces the inputs state.
export function setInputs(inputs: SessionStatusResponse['session']['inputs']) {
    return actionObject('SET_INPUTS', {inputs});
}

// Adds any inputs specified, but does not fully replace inputs.
// Removes the inputs by the given source name from the inputs and session states
export function removeInputs(inputs: string[]) {
    return actionObject('REMOVE_INPUTS', {inputs});
}

// Removes the downloads by the given step names from the status.  Usually this happens from server
// responses anyways, but it is a responsive UX feature.
export function removeDownloads(stepNames: string[]) {
    return actionObject('REMOVE_DOWNLOADS', {stepNames});
}

export function addValidationErrors(inputName: string, inputLabel: string, errors: string[]) {
    return actionObject('ADD_VALIDATION_ERRORS', {inputName, inputLabel, errors});
}

export function removeValidationErrors(inputName: string) {
    return actionObject('REMOVE_VALIDATION_ERRORS', {inputName});
}

export function startPolling() {
    return actionObject('MODIFY_POLLING', {delta: 1})
}

export function finishPolling() {
    return actionObject('MODIFY_POLLING', {delta: -1})
}

// value = null => no value is set, value = [null] => the value null has been set.
export function setBufferedInput(source: string, value: null | [any]) {
    return actionObject('SET_BUFFERED_INPUT', {source, value});
}

export const clearBufferedInput = actionObject('CLEAR_BUFFERED_INPUT', {});

export type VulcanAction = ReturnType<typeof setWorkflows> | ReturnType<typeof setWorkflow> | ReturnType<typeof setStatus> |
    ReturnType<typeof setDownloadedData> | ReturnType<typeof setSession> | ReturnType<typeof setInputs> |
    ReturnType<typeof releaseDownloadedData> | ReturnType<typeof removeInputs> |
    ReturnType<typeof removeDownloads> |
    ReturnType<typeof startPolling> | ReturnType<typeof finishPolling> |
    ReturnType<typeof addValidationErrors> | ReturnType<typeof removeValidationErrors> |
    ReturnType<typeof setBufferedInput> | typeof clearBufferedInput;
