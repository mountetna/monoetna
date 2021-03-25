// Status
export const STATUS = {
  PENDING: 'pending',
  COMPLETE: 'complete',
  ERROR: 'error',
  RUNNING: 'running',
};

// Types
export const TYPE = {
  INTEGER: 'int',
  FLOAT: 'float',
  BOOL: 'boolean',
  STRING: 'string',
  ARRAY: 'array',
  FILE: 'File',
  MULTISELECT_STRING: 'multiselect-string',
  SELECT_AUTOCOMPLETE: 'select-autocomplete',
  CHECKBOXES: 'checkboxes',
  NESTED_SELECT_AUTOCOMPLETE: 'nested-select-autocomplete',
};

// CWL Step RUN Sentinels
export const RUN = {
  UI_QUERY: 'ui-queries/',
  UI_OUTPUT: 'ui-outputs/'
};

// UI Output widgets
export const OUTPUT_COMPONENT = {
  LINK: 'link',
  PLOTLY: 'plotly',
  CONSIGNMENT: 'consignment',
  RAW: 'raw'
};

export const defaultWorkflowsResponse = {
  workflows: [] as Workflow[],
}

export type WorkflowsResponse = typeof defaultWorkflowsResponse;

export interface StepInput {
  id: string,
  source: string,
  doc?: string,
  label?: string,
}

export interface WorkflowStep {
  name: string,
  run: string,
  in: StepInput[],
  out: string[],
  label: string | undefined,
  doc: string | undefined,
}

export interface WorkflowOutput {
  outputSource: string,
  label: string | undefined,
  type: string,
  default: any | null,
  format: string | null,
}

export interface WorkflowInput {
  doc: string | undefined,
  label: string | undefined,
  type: string,
  format: string | null,
  default: any | null,
}

export interface Workflow {
  name: string,
  inputs: {[k: string]: WorkflowInput},
  outputs: {[k: string]: WorkflowOutput},
  dependencies_of_outputs: {[k: string]: string[]},
  steps: [WorkflowStep[]],
  vignette?: string,
  image?: string,
  projects?: string[],
  authors?: string[],
}

export type StatusString = 'running' | 'pending' | 'complete' | 'error';

export interface StepStatus {
  name: string,
  status: StatusString,
  downloads: {[k: string]: string} | undefined | null,
  error: string | undefined | null,
}

export const defaultVulcanSession = {
  project_name: "",
  workflow_name: "",
  key: "",
  inputs: {} as {[k: string]: any},
};

export type VulcanSession = typeof defaultVulcanSession;


export const defaultSessionStatusResponse = {
  outputs: {} as  {[k: string]: {downloads: StepStatus['downloads'], status: StepStatus['status']}},
  session: defaultVulcanSession,
  status: [[]] as [StepStatus[]]
}

export type SessionStatusResponse = typeof defaultSessionStatusResponse;
