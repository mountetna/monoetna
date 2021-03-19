import {Workflow, WorkflowsResponse} from "../api/types";

export const SET_WORKFLOWS = 'SET_WORKFLOWS';
export const SET_WORKFLOW = 'SET_WORKFLOW';
export const SET_STATUS = 'SET_STATUS';
export const SET_DATA = 'SET_DATA';
export const SET_PATH = 'SET_PATH';
export const SET_STEP = 'SET_STEP';
export const SET_SESSION = 'SET_SESSION';
export const SET_INPUTS = 'SET_INPUTS';
export const COMMIT_INPUTS = 'COMMIT_INPUTS';

function actionObject<T extends string, P>(type: T, payload: P): { type: T } & P {
    return { ...payload, type };
}

export function setWorkflows(workflows: WorkflowsResponse['workflows']) {
    return actionObject('SET_WORKFLOWS', {workflows});
}

export function setWorkflow(workflow: Workflow) {
    return actionObject('SET_WORKFLOW', {workflow});
}

export type VulcanAction = ReturnType<typeof setWorkflows> | ReturnType<typeof setWorkflow>;
