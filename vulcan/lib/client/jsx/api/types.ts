export type Label = string | null;

export interface WorkflowsResponse {
    workflows: Workflow[]
}

export interface StepInput {
    id: string,
    source: [string, string],
}

export interface WorkflowStep {
    name: string,
    run: string,
    in: StepInput[],
    out: string[],
    label: Label,
    doc: string | undefined,
}

export interface WorkflowOutput {
    outputSource: string,
    label: Label,
    type: string,
    default: any | null,
    format: string | null,
}

export interface WorkflowInput {
    label: Label,
    type: string,
    format: string | null,
    default: any | null,
}

export interface Workflow {
    name: string,
    inputs: {[k: string]: WorkflowInput},
    outputs: {[k: string]: WorkflowOutput},
    steps: [WorkflowStep[]],
    vignette?: string,
    image?: string,
    projects?: string[],
    authors?: string[],
}

export type StatusString = 'running' | 'pending' | 'complete' | 'error'
export interface StepStatus {
    name: string,
    status: StatusString,
    downloads: {[k: string]: string} | undefined | null,
    error: string | undefined | null,
}

export interface SessionStatusResponse {
    outputs: {[k: string]: {downloads: StepStatus['downloads'], status: StepStatus['status']}},
    session: VulcanSession,
    status: [StepStatus[]]
}

export interface VulcanSession {
    project_name: string,
    workflow_name: string,
    key: string,
    inputs: {[k: string]: any}
}


