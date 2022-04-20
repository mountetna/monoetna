import {
  SessionStatusResponse,
  VulcanFigure,
  VulcanFigureSession,
  VulcanSession,
  Workflow,
  WorkflowsResponse
} from '../api_types';
import {Maybe} from '../selectors/maybe';
import {selectFigure, selectSession} from '../selectors/workflow_selectors';

function actionObject<T extends string, P>(type: T, payload: P): {type: T} & P {
  return {...payload, type};
}

export function setWorkflows(workflows: WorkflowsResponse['workflows']) {
  return actionObject('SET_WORKFLOWS', {workflows});
}

export function setWorkflow(workflow: Workflow, projectName: string) {
  return actionObject('SET_WORKFLOW', {workflow, projectName});
}

export function setStatus(
  status: SessionStatusResponse['status'],
  submittingStep: Maybe<string> = null
) {
  return actionObject('SET_STATUS', {status, submittingStep});
}

export function setDownloadedData(url: string, data: any) {
  return actionObject('SET_DOWNLOAD', {url, data});
}

export function setSession(session: SessionStatusResponse['session']) {
  return actionObject('SET_SESSION', {session});
}

export function setSessionAndFigure(figureSession: VulcanFigureSession) {
  return setSessionAndFigureSeparately(
    selectFigure(figureSession),
    selectSession(figureSession)
  );
}

export function setSessionAndFigureSeparately(
  figure: VulcanFigure,
  session: VulcanSession
) {
  return actionObject('SET_SESSION_AND_FIGURE', {figure, session});
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

export function setStepUIAutoPass(step: string | null) {
  return actionObject('SET_AUTO_PASS_STEP', {step});
}

export function clearStepUIAutoPass(step: string | null) {
  return actionObject('CLEAR_AUTO_PASS_STEP', {step});
}

export function clearCommitTrigger(step: string | null) {
  return actionObject('CLEAR_COMMIT_TRIGGER', {step});
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
  | ReturnType<typeof setWorkflows>
  | ReturnType<typeof setWorkflow>
  | ReturnType<typeof setStatus>
  | ReturnType<typeof setDownloadedData>
  | ReturnType<typeof setSession>
  | ReturnType<typeof setInputs>
  | ReturnType<typeof removeInputs>
  | ReturnType<typeof removeDownloads>
  | ReturnType<typeof startPolling>
  | ReturnType<typeof finishPolling>
  | ReturnType<typeof addValidationErrors>
  | ReturnType<typeof removeValidationErrors>
  | ReturnType<typeof setBufferedInput>
  | ReturnType<typeof clearBufferedInput>
  | ReturnType<typeof setStepUIAutoPass>
  | ReturnType<typeof clearStepUIAutoPass>
  | ReturnType<typeof clearCommitTrigger>
  | ReturnType<typeof setRunTrigger>
  | ReturnType<typeof clearRunTriggers>
  | ReturnType<typeof checkCommittedStepPending>
  | ReturnType<typeof clearCommittedStepPending>
  | ReturnType<typeof setSessionAndFigure>
  | ReturnType<typeof setSessionAndFigureSeparately>;
