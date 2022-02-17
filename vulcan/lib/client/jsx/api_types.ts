// Status
export const STATUS = {
  PENDING: 'pending',
  COMPLETE: 'complete',
  ERROR: 'error',
  RUNNING: 'running'
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
  MULTIPLE_MULTISELECT_STRING_ALL: 'multiple-multiselect-string-all',
  MULTIPLE_STRING: 'multiple-string',
  SINGLE_DROPDOWN_MULTICHECKBOX: 'single-dropdown-multicheckbox',
  SCATTER_PLOTLY: 'scatter-plotly',
  BAR_PLOTLY: 'bar-plotly',
  Y_PLOTLY: 'y-plotly',
  DIFF_EXP_SC: 'diff-exp-sc'
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
  workflows: [] as Workflow[]
};

export type WorkflowsResponse = typeof defaultWorkflowsResponse;
export type ProjectWorkflowsResponse = typeof defaultWorkflowsResponse;

export interface StepInput {
  id: string;
  source: string;
  doc?: string | null;
  label?: string | null;
}

export interface WorkflowStep {
  name: string;
  run: string;
  in: StepInput[];
  out: string[];
  label?: string | null;
  doc?: string | null;
}

export const defaultWorkflowStep: WorkflowStep = {
  name: '',
  run: '',
  in: [],
  out: []
};

export interface WorkflowOutput {
  outputSource: string;
  label?: string | null;
  type: string;
  default: any | null;
  format: string | null;
}

export const defaultWorkflowInput: WorkflowInput = {
  type: TYPE.STRING
};

export interface WorkflowInput {
  doc?: string | null;
  label?: string | null;
  type: string;
  format?: string | null;
  default?: any | null;
}

export interface Workflow {
  class: string;
  cwlVersion: string;
  name: string;
  inputs: {[k: string]: WorkflowInput};
  outputs: {[k: string]: WorkflowOutput};
  steps: [WorkflowStep[]];
  vignette?: string;
  image?: string;
  projects?: string[];
  authors?: string[];
  displayName?: string;
  description?: string;
  tags?: string[];
}

export const defaultWorkflow: Workflow = {
  class: '',
  cwlVersion: '',
  name: '',
  inputs: {} as {[k: string]: WorkflowInput},
  outputs: {} as {[k: string]: WorkflowOutput},
  steps: [[]] as [WorkflowStep[]]
};

export type StatusString = 'running' | 'pending' | 'complete' | 'error';

export interface StepStatus {
  name: string;
  status: StatusString;
  downloads?: {[k: string]: string} | null;
  error?: string | null;
  hash: string;
}

export const defaultStepStatus: StepStatus = {
  name: '',
  status: 'pending',
  hash: ''
};

export const defaultVulcanSession = {
  project_name: '',
  workflow_name: '',
  key: '',
  inputs: {} as {[k: string]: any}
};

export type VulcanSession = typeof defaultVulcanSession;

export const defaultSessionStatusResponse = {
  outputs: {downloads: null, status: 'pending'} as {
    downloads: StepStatus['downloads'];
    status: StepStatus['status'];
  },
  session: defaultVulcanSession,
  status: [[]] as [StepStatus[]]
};

export type SessionStatusResponse = typeof defaultSessionStatusResponse;

export interface VulcanFigure {
  figure_id: number | null;
  title?: string | null;
  author?: string;
}

export type VulcanFigureResponse = VulcanSession & VulcanFigure;

export const defaultFigure = {
  figure_id: null
}