import {setStatus, setWorkflow, setWorkflows, VulcanAction} from "../actions/vulcan_actions";
import {workflowsResponse} from "./fixtures/workflows-response";
import {statusWithoutDownloads} from "./fixtures/status-without-downloads";
import {
    defaultFigure,
    defaultSessionStatusResponse,
    defaultStepStatus,
    defaultWorkflowStep,
    SessionStatusResponse,
    StepStatus, VulcanFigure, VulcanSession,
    Workflow,
    WorkflowsResponse,
    WorkflowStep
} from "../api_types";
import {defaultSession, VulcanState} from "../reducers/vulcan_reducer";

export function createWorkflowFixture(workflow: Partial<Workflow>): Workflow {
    return { ...workflowsResponse['workflows'][0], projects: ['test'], ...workflow };
}

export function createStepFixture(step: Partial<WorkflowStep>): WorkflowStep {
    return { ...defaultWorkflowStep, ...step };
}

export function createStepStatusFixture(stepStatus: Partial<StepStatus>): StepStatus {
    return {...defaultStepStatus, ...stepStatus};
}

export function createStatusFixture(workflow: Workflow, ...statusPatches: Partial<StepStatus>[]): [StepStatus[]] {
    return [workflow.steps[0].map(step => {
      const patch = statusPatches.find(({name}) => name === step.name);
      if (patch) return {...defaultStepStatus, ...patch};
      return {...defaultStepStatus, name: step.name };
    })];
}

export function createSessionFixture(workflow_name: string, session: Partial<VulcanSession> = {}): VulcanSession {
    return {
        ...defaultSession,
        ...session,
        workflow_name,
    }
}

export function createStoredSessionFixture(workflow_name: string, session: Partial<VulcanSession> = {}, figure: Partial<VulcanFigure> = {}): VulcanSession & VulcanFigure {
    return {
        ...defaultSession,
        ...defaultFigure,
        ...session,
        ...figure,
        workflow_name,
    }
}

export function createUpdatedStatusFixture(workflow: Workflow, status: VulcanState['status'], ...statusPatches: Partial<StepStatus>[]): [StepStatus[]] {
    return [workflow.steps[0].map(step => {
        const patch = statusPatches.find(({name}) => name === step.name);
        if (patch) return {...defaultStepStatus, ...patch};
        const original = status[0].find(({name}) => name === step.name);
        if (original) return original;
        return {...defaultStepStatus, name: step.name };
    })];
}

export function createStatusResponseFixture(response: Partial<SessionStatusResponse>): SessionStatusResponse {
    return {...defaultSessionStatusResponse, ...response};
}

export function createUpdatedStatusResponseFixture(state: VulcanState, response: Partial<SessionStatusResponse>): SessionStatusResponse {
    return {
        ...defaultSessionStatusResponse,
        session: {...state.session, ...(response.session || {})},
        status: response.status || state.status,
    };
}

export function findWorkflowFromResponse(response: WorkflowsResponse, n: string): Workflow {
    const result = response.workflows.find(({name}) => name === n);
    if (!result) throw new Error("Workflow named " + n + " not found on the response");
    return result;
}