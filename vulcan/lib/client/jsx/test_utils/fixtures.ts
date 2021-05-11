import {setStatus, setWorkflow, setWorkflows, VulcanAction} from "../actions/vulcan";
import {workflowsResponse} from "./fixtures/workflows-response";
import {statusWithoutDownloads} from "./fixtures/status-without-downloads";
import {
    defaultSessionStatusResponse,
    defaultStepStatus,
    defaultWorkflowStep,
    SessionStatusResponse,
    StepStatus,
    Workflow, WorkflowsResponse,
    WorkflowStep
} from "../api_types";

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

export function createStatusResponseFixture(response: Partial<SessionStatusResponse>): SessionStatusResponse {
    return {...defaultSessionStatusResponse, ...response};
}

export function findWorkflowFromResponse(response: WorkflowsResponse, n: string): Workflow {
    const result = response.workflows.find(({name}) => name === n);
    if (!result) throw new Error("Workflow named " + n + " not found on the response");
    return result;
}