import {
  AccountingReturn,
  SessionStatusResponse,
  VulcanFigure,
  VulcanFigureSession,
  VulcanSession,
  Workflow,
  WorkflowsResponse,
  Workspaces,
  Workspace,
  WorkspaceStatus
} from '../api_types';
import {Maybe} from '../selectors/maybe';
import { selectFigure, selectSession, workspaceId } from '../selectors/workflow_selectors';

function actionObject<T extends string, P>(type: T, payload: P): {type: T} & P {
  return {...payload, type};
}

export function setProject(project: string) {
  return actionObject('SET_PROJECT', {project});
}

export function setWorkflows(workflows: WorkflowsResponse) {
  return actionObject('SET_WORKFLOWS', {workflows});
}

export function setWorkspaces(workspaces: Workspaces) {
  return actionObject('SET_WORKSPACES', {workspaces});
}

export function setWorkflow(workflow: Workflow, projectName: string) {
  return actionObject('SET_WORKFLOW', {workflow, projectName});
}

export function setWorkspaceId(workspaceId: Workspace['workspace_id']) {
  return actionObject('SET_WORKSPACE_ID', {workspaceId});
}

export function setWorkspace(workspace: Workspace, projectName: string) {
  return actionObject('SET_WORKSPACE', {workspace, projectName});
}

export function setConfigId(configId: Workspace['workspace_id']) {
  return actionObject('SET_CONFIG_ID', {configId});
}

export function setRunId(runId: Workspace['workspace_id']) {
  return actionObject('SET_RUN_ID', {runId});
}

export function setStatus(
  accounting: AccountingReturn,
  submittingStep: Maybe<string> = null
) {
  return actionObject('SET_STATUS', {accounting, submittingStep});
}

export function setDownloadedData(fileName: string, fileData: any) {
  return actionObject('SET_DOWNLOAD', {fileName, fileData});
}

// Fully replaces the inputs state.
export function setUIValues(values: WorkspaceStatus['ui_contents']) {
  return actionObject('SET_UI_VALUES', {values});
}

// // Removes the downloads by the given step names from the status.  Usually this happens from server
// // responses anyways, but it is a responsive UX feature.
// export function removeDownloads(stepNames: string[]) {
//   return actionObject('REMOVE_DOWNLOADS', {stepNames});
// }

export function addValidationErrors(
  stepName: string | null,
  inputLabel: string,
  errors: string[]
) {
  return actionObject('ADD_VALIDATION_ERRORS', {stepName, inputLabel, errors});
}

// Relies on the object identity of the errors string, but loosens the necessary uniqueness of inputLabel
export function removeValidationErrors(errors: string[]) {
  return actionObject('REMOVE_VALIDATION_ERRORS', {errors});
}

export function startPolling() {
  return actionObject('MODIFY_POLLING', {delta: 1});
}

export function finishPolling() {
  return actionObject('MODIFY_POLLING', {delta: -1});
}

export function setBufferedInput(step: string | null) {
  return actionObject('SET_BUFFERED_INPUT', {step});
}

export function clearBufferedInput(step: string | null) {
  return actionObject('CLEAR_BUFFERED_INPUT', {step});
}

export function setAutoPassStep(step: string | null) {
  return actionObject('SET_AUTO_PASS_STEP', {step});
}

export function clearAutoPassStep(step: string | null) {
  return actionObject('CLEAR_AUTO_PASS_STEP', {step});
}

export function setRunTrigger(step: string | null) {
  return actionObject('SET_RUN_TRIGGER', {step});
}

export function clearRunTriggers(steps: (string | null)[]) {
  return actionObject('CLEAR_RUN_TRIGGERS', {steps});
}

export function checkCommittedStepPending(step: Maybe<string>) {
  return actionObject('CHECK_CHANGES_READY', {step});
}

export function clearCommittedStepPending() {
  return actionObject('CLEAR_CHANGES_READY', {});
}

export type VulcanAction =
  | ReturnType<typeof setProject>
  | ReturnType<typeof setWorkflows>
  | ReturnType<typeof setWorkspaces>
  | ReturnType<typeof setWorkflow>
  | ReturnType<typeof setWorkspace>
  | ReturnType<typeof setWorkspaceId>
  | ReturnType<typeof setConfigId>
  | ReturnType<typeof setRunId>
  | ReturnType<typeof setStatus>
  | ReturnType<typeof setDownloadedData>
  | ReturnType<typeof setUIValues>
  // | ReturnType<typeof removeDownloads>
  | ReturnType<typeof startPolling>
  | ReturnType<typeof finishPolling>
  | ReturnType<typeof addValidationErrors>
  | ReturnType<typeof removeValidationErrors>
  | ReturnType<typeof setBufferedInput>
  | ReturnType<typeof clearBufferedInput>
  | ReturnType<typeof setAutoPassStep>
  | ReturnType<typeof clearAutoPassStep>
  | ReturnType<typeof setRunTrigger>
  | ReturnType<typeof clearRunTriggers>
  | ReturnType<typeof checkCommittedStepPending>
  | ReturnType<typeof clearCommittedStepPending>
