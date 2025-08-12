import {
  AccountingReturn,
  VulcanStorage,
  Workflow,
  Workspace,
  WorkspaceStatus,
  WorkspaceRaw,
  RunStatus,
  MultiFileContent,
} from '../api_types';
import { DataEnvelope } from '../components/workspace/ui_definitions/input_types';
import {Maybe} from '../selectors/maybe';

function actionObject<T extends string, P>(type: T, payload: P): {type: T} & P {
  return {...payload, type};
}

export function setWorkflow(workflow: Workflow, projectName: string) {
  return actionObject('SET_WORKFLOW', {workflow, projectName});
}

export function setWorkflowsWorkspaces(updates: any) {
  return actionObject('SET_WORKFLOWS_AND_WORKSPACES', {updates});
}

export function updateWorkflowsWorkspaces() {
  return actionObject('UPDATE_WORKFLOWS_AND_WORKSPACES', {});
}

export function setWorkspaceId(workspaceId: Workspace['workspace_id']) {
  return actionObject('SET_WORKSPACE_ID', {workspaceId});
}

export function setWorkspace(workspace: Workspace | WorkspaceRaw, projectName: string) {
  return actionObject('SET_WORKSPACE', {workspace, projectName});
}

export function setFullWorkspaceState(workspace: Workspace, status: WorkspaceStatus, update_files: boolean, isRunning: boolean ) {
  return actionObject('SET_FULL_WORKSPACE_STATUS', {workspace, status, update_files, isRunning});
}

export function setStateFromStorage(storage: VulcanStorage) {
  return actionObject('SET_STATE_FROM_STORAGE', {storage});
}

export function setWorkspaceStateSyncd() {
  return actionObject('SET_WORKSPACE_STATE_SYNCD', {});
}

export function setConfigId(configId: Workspace['workspace_id']) {
  return actionObject('SET_CONFIG_ID', {configId});
}
export function setAttemptingToRun(to: boolean) {
  return actionObject('SET_ATTEMPTING_RUN', {to});
}

export function setAttemptingToCancel(to: boolean) {
  return actionObject('SET_ATTEMPTING_CANCELLATION', {to});
}

export function setRunning(runId: Workspace['workspace_id']) {
  return actionObject('SET_RUNNING', {runId});
}

export function setLastConfig(lastConfig: WorkspaceStatus['last_params']) {
  return actionObject('SET_LAST_CONFIG', {lastConfig})
}

export function useUIAccounting(
  accounting: AccountingReturn,
  submittingStep: Maybe<string> = null,
) {
  return actionObject('USE_UI_ACCOUNTING', {accounting, submittingStep});
}

export function setStatusFromStatuses(
  statusReturns: RunStatus,
  isRunning: boolean
) {
  return actionObject('SET_STATUS_FROM_STATUSES', {statusReturns, isRunning});
}

export function setWorkspaceFiles(fileNames: string[]) {
  return actionObject('SET_WORKSPACE_FILES', {fileNames});
}

// Fully replaces the inputs state.
export function setUIValues(values: DataEnvelope<Maybe<any>>, stepName: string) {
  return actionObject('SET_UI_VALUES', {values, stepName});
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

export function startSyncing() {
  return actionObject('MODIFY_SYNCING', {to: true, submittingStep: undefined});
}

export function endSyncing(submittingStep: string) {
  return actionObject('MODIFY_SYNCING', {to: false, submittingStep: submittingStep});
}

export function updateFiles(statusUpdates: Pick<WorkspaceStatus, 'output_files' | 'file_contents'>) {
  return actionObject('UPDATE_FILES', {statusUpdates});
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

export function checkCommittedStepPending(step: string | null) {
  return actionObject('CHECK_CHANGES_READY', {step});
}

export function clearRunning() {
  return actionObject('CLEAR_RUNNING', {});
}

export function clearCommittedStepPending() {
  return actionObject('CLEAR_CHANGES_READY', {});
}

export type VulcanAction =
  | ReturnType<typeof setWorkflow>
  | ReturnType<typeof setWorkflowsWorkspaces>
  | ReturnType<typeof updateWorkflowsWorkspaces>
  | ReturnType<typeof setWorkspace>
  | ReturnType<typeof setWorkspaceId>
  | ReturnType<typeof setFullWorkspaceState>
  | ReturnType<typeof setStateFromStorage>
  | ReturnType<typeof setWorkspaceStateSyncd>
  | ReturnType<typeof setConfigId>
  | ReturnType<typeof setAttemptingToRun>
  | ReturnType<typeof setRunning>
  | ReturnType<typeof setAttemptingToCancel>
  | ReturnType<typeof setLastConfig>
  | ReturnType<typeof useUIAccounting>
  | ReturnType<typeof setStatusFromStatuses>
  | ReturnType<typeof setWorkspaceFiles>
  | ReturnType<typeof setUIValues>
  | ReturnType<typeof startSyncing>
  | ReturnType<typeof endSyncing>
  | ReturnType<typeof updateFiles>
  | ReturnType<typeof addValidationErrors>
  | ReturnType<typeof removeValidationErrors>
  | ReturnType<typeof setBufferedInput>
  | ReturnType<typeof clearBufferedInput>
  | ReturnType<typeof setAutoPassStep>
  | ReturnType<typeof clearAutoPassStep>
  | ReturnType<typeof setRunTrigger>
  | ReturnType<typeof clearRunTriggers>
  | ReturnType<typeof clearRunning>
  | ReturnType<typeof checkCommittedStepPending>
  | ReturnType<typeof clearCommittedStepPending>
