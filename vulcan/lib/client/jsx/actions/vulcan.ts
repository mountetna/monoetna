import {SessionStatusResponse, Workflow, WorkflowsResponse} from "../api_types";

function actionObject<T extends string, P>(type: T, payload: P): { type: T } & P {
    return { ...payload, type };
}

export function setWorkflows(workflows: WorkflowsResponse['workflows']) {
    return actionObject('SET_WORKFLOWS', {workflows});
}

export function setWorkflow(workflow: Workflow) {
    return actionObject('SET_WORKFLOW', {workflow});
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

export function setInputs(inputs: SessionStatusResponse['session']['inputs']) {
    return actionObject('SET_INPUTS', {inputs});
}

export function commitInputs() {
    return actionObject('COMMIT_INPUTS', {});
}

export function setCalculating(calculating: boolean) {
    return actionObject('SET_CALCULATING', {calculating});
}

const actions = {
    commitInputs,
}

export type VulcanAction = ReturnType<typeof setWorkflows> | ReturnType<typeof setWorkflow> | ReturnType<typeof setStatus> |
    ReturnType<typeof setDownloadedData> | ReturnType<typeof setSession> | ReturnType<typeof setInputs> |
    ReturnType<typeof commitInputs> | ReturnType<typeof releaseDownloadedData> | ReturnType<typeof setCalculating>;
