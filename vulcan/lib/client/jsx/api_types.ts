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
  METIS_FILE: 'MetisFile',
  METIS_CSV_OR_TSV: 'MetisCSVorTSV',
  METIS_FOLDER: 'MetisFolder',
  METIS_FILE_OR_FOLDER: 'MetisPath',
  MULTISELECT_STRING: 'multiselect-string',
  SELECT_AUTOCOMPLETE: 'select-autocomplete',
  SELECT_AUTOCOMPLETE_MULTI_PICK: 'select-autocomplete-multi-pick',
  CHECKBOXES: 'checkboxes',
  NESTED_SELECT_AUTOCOMPLETE: 'nested-select-autocomplete',
  NESTED_SELECT_AUTOCOMPLETE_MULTI_PICK: 'nested-select-autocomplete-multi-pick',
  MULTIPLE_MULTISELECT_STRING_ALL: 'multiple-multiselect-string-all',
  MULTIPLE_STRING: 'multiple-string',
  SINGLE_DROPDOWN_MULTICHECKBOX: 'single-dropdown-multicheckbox',
  SCATTER_PLOTLY: 'scatter-plotly',
  BAR_PLOTLY: 'bar-plotly',
  Y_PLOTLY: 'y-plotly',
  DITTOSEQ_DIM_PLOT: 'dittoseq-dim-plot',
  DITTOSEQ_SCATTER_PLOT: 'dittoseq-scatter-plot',
  DITTOSEQ_BAR_PLOT: 'dittoseq-bar-plot',
  DITTOSEQ_PLOT: 'dittoseq-plot',
  ANY_DITTOSEQ: 'any-dittoseq',
  DIFF_EXP_SC: 'diff-exp-sc',
  DATA_TRANSFORMATION: 'data-transformation',
  ANNOTATION_EDITOR: 'annotation-editor',
  ANY_VIZ: 'any-viz'
};

// CWL Step RUN Sentinels
export const RUN = {
  UI_QUERY: 'ui-queries/',
  UI_OUTPUT: 'ui-outputs/'
};

// UI Output widgets
export const OUTPUT_COMPONENT = {
  LINK: 'link',
  PLOT: 'plot',
  PLOTLY: 'plotly',
  PNG: 'png',
  CONSIGNMENT: 'consignment',
  RAW: 'raw'
};

export const defaultWorkflowsResponse = [] as Workflow[];

export type WorkflowsResponse = typeof defaultWorkflowsResponse;
// export type ProjectWorkflowsResponse = typeof defaultWorkflowsResponse;

// export interface StepInput {
//   id: string;
//   source: string;
//   doc?: string | null;
//   label?: string | null;
// }

export interface WorkspaceStep {
  name: string;
  // run: string;
  vulcan_config: boolean;
  in: StepInput[];
  out: string[];
  label?: string | null;
  doc?: string | null;
}

// export const defaultWorkflowStep: WorkflowStep = {
//   name: '',
//   run: '',
//   in: [],
//   out: []
// };

// export interface WorkflowOutput {
//   outputSource: string;
//   label?: string | null;
//   type: string;
//   default: any | null;
//   format: string | null;
// }

// export const defaultWorkflowInput: WorkflowInput = {
//   type: TYPE.STRING
// };

// export interface WorkflowInput {
//   doc?: string | null;
//   label?: string | null;
//   type: string;
//   format?: string | null;
//   default?: any | null;
// }

export interface Workflow {
  id: number;
  name: string;
  projects: string[];
  branch: string;
  repo_remote_url: string;
  created_at: number;
  updated_at: number;
  vignette?: string;
  image?: string;
  authors?: string[];
  displayName?: string;
  description?: string;
  tags?: string[];
  lastModified?: string;
}

export const defaultWorkflow: Workflow = {
  id: 0,
  name: '',
  projects: [],
  branch: '',
  repo_remote_url: '',
  created_at: 0,
  updated_at: 0
};

export interface Workspace {
  workspace_id: number;
  workflow_id: number;
  project: string;
  vulcan_config: VulcanConfig;
  dag: string[];
  last_config?: {[k: string]: any}
  last_job_status?: {[k: string]: StatusString}
}

export type Workspaces = Workspace[]

export interface FileContentResponse {
  filename: string;
  content: string;
}

export type MultiFileContentResponse = FileContentResponse[]

export interface MultiFileContent {
  [filename: string]: any;
}

export interface ParamsContent {
  [k: string]: any
}

export interface AccountingReturn {
  config_id: number,
  scheduled: string[],
  downstream: string[],
}

// export type VulcanConfig = VulcanConfigElement[]
export type VulcanConfig = {[k: string]: VulcanConfigElement}

export interface VulcanConfigElement {
  display: string;
  ui_component: string;
  input?: InputConfig;
  output?: OutputConfig;
}

export interface InputConfig {
  files?: string[]
  params?: string[]
}

export interface OutputConfig {
  files?: string[]
  params?: string[]
}

export interface RunReturn {
  run_id: number
}

export interface RunStatus {
  [k: string]: StatusString
}

export type sacctStatusString =
  'BOOT_FAIL' |
  'CANCELLED' |
  'COMPLETED' |
  'CONFIGURING' |
  'COMPLETING' |
  'DEADLINE' |
  'FAILED' |
  'NODE_FAIL' |
  'OUT_OF_MEMORY' |
  'PENDING' |
  'PREEMPTED' |
  'REBOOT_PENDING' |
  'REQUEUED' |
  'RESIZING' |
  'REVOKED' |
  'RUNNING' |
  'SUSPENDED' |
  'TIMEOUT'

export type StatusStringFine = sacctStatusString | 'NOT STARTED'

export const StatusStringBroaden = (fine: StatusStringFine) => {
  return ({
    'BOOT_FAIL': 'error',
    'CANCELLED': 'error',
    'COMPLETED': 'complete',
    'CONFIGURING': 'running',
    'COMPLETING': 'running',
    'DEADLINE': 'error',
    'FAILED': 'error',
    'NODE_FAIL': 'error',
    'OUT_OF_MEMORY': 'error',
    'PENDING': 'upcoming',
    'PREEMPTED': 'error',
    'REBOOT_PENDING': 'upcoming',
    'REQUEUED': 'upcoming',
    'RESIZING': 'running',
    'REVOKED': 'error',
    'RUNNING': 'running',
    'SUSPENDED': 'error',
    'TIMEOUT': 'error',
    'NOT STARTED': 'pending'
  } as {[k:string]: StatusString})[fine]
}

export type StatusString = 'running' | 'upcoming' |  'pending' | 'complete' | 'error';

export interface StepStatus {
  name: string;
  status: StatusString;
  statusFine: StatusStringFine;
  outputs?: {files?: string[], params?: string[]};
  error?: string;
}

export const defaultStepStatus: StepStatus = {
  name: '',
  status: 'pending',
  statusFine: 'NOT STARTED'
};

export type WorkspaceStatus = {
  steps: {[k: string]: StepStatus},
  output_files: string[],
  config_contents: {[k: string]: any},
  ui_contents: {[k: string]: {[k: string]: any}},
}

export const defaultWorkspaceStatus: WorkspaceStatus = {
  steps: {},
  config_contents: {},
  ui_contents: {},
  output_files: []
}

export const defaultVulcanSession = {
  project_name: '',
  workflow_id: null as number | null,
  workspace_id: null as number | null,
  workspace: null as Workspace | null, 
  key: '',
  params: {} as {[k: string]: any},
  ui_values: {} as {[k: string]: any},
};

export type VulcanSession = typeof defaultVulcanSession;

// export const defaultSessionStatusResponse = {
//   session: defaultVulcanSession,
//   status: [[]] as [StepStatus[]],
//   files: [] as string[],
// };

// export type SessionStatusResponse = typeof defaultSessionStatusResponse;

// // Update me!
// export interface VulcanFigure {
//   id: number | null;
//   figure_id?: number | null;
//   inputs: {[k: string]: any};
//   title?: string;
//   author?: string;
//   thumbnails?: string[];
//   comment?: string;
//   tags?: string[];
//   workflow_snapshot?: Workflow;
// }

// export type VulcanFigureSession = VulcanSession & VulcanFigure;

// export interface VulcanRevision {
//   inputs: {[k: string]: any};
//   title?: string;
//   tags?: string[];
//   id: number;
//   workflow_snapshot?: Workflow;
//   dependencies: {[key: string]: string};
// }

// export const defaultFigure = {
//   id: null,
//   figure_id: null,
//   inputs: {},
//   workflow_snapshot: defaultWorkflow
// };

// export interface FiguresResponse {
//   figures: VulcanFigureSession[];
// }
