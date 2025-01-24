import { Maybe } from "./selectors/maybe";

// CWL Step RUN Sentinels
export const RUN = {
  UI_QUERY: 'ui-queries/',
  UI_OUTPUT: 'ui-outputs/'
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
  input: InputOutputConfig;
  output: InputOutputConfig;
  vulcan_config?: boolean;
  label?: string;
  ui_component?: string;
  doc?: string;
}

export const defaultWorkspaceStep: WorkspaceStep = {
  name: '',
  input: {},
  output: {}
};

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
  id: number | null;
  name: string;
  project: string;
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
  id: null,
  name: '',
  project: '',
  branch: '',
  repo_remote_url: '',
  created_at: 0,
  updated_at: 0,
  tags: []
};

interface WorkspaceMinusInconsistent {
  workspace_id: number | null;
  workflow_id: number | null;
  project?: string;
  steps: {[k: string]: WorkspaceStep};
  dag: string[];
  last_config?: {[k: string]: any};
  last_job_status?: {[k: string]: StatusStringFine};
  vignette?: string;
  thumbnails?: string[];
  author?: string;
  name?: string;
  tags: string[];
}

export interface Workspace extends WorkspaceMinusInconsistent {
  vulcan_config: VulcanConfig;
}

export interface WorkspaceRaw extends WorkspaceMinusInconsistent {
  vulcan_config: VulcanConfigElement[]
}

export const defaultWorkspace: Workspace = {
  workspace_id: null,
  workflow_id: null,
  vulcan_config: {},
  steps: {},
  dag: [],
  last_config: {},
  last_job_status: {},
  tags: [],
};

export type Workspaces = Workspace[]

export type WorkspacesResponse = {workspaces: Workspaces}

export interface FileContentResponse {
  filename: string;
  content: string;
}

export type MultiFileContentResponse = FileContentResponse[]

export interface MultiFileContent {
  [filename: string]: any;
}

export type FlatParams = {
  [k: string]: any
};

export interface AccountingReturn {
  config_id: number,
  scheduled: string[],
  downstream: string[],
}

// export type VulcanConfig = VulcanConfigElement[]
export type VulcanConfig = {[k: string]: VulcanConfigElement}
export type VulcanConfigRaw = VulcanConfigElement[];

export type VulcanConfigElement = {
  name: string;
  display: string;
  ui_component: string;
  default?: any;
  doc?: string;
  input?: InputOutputConfig;
  input_map?: string[];
  output?: InputOutputConfig;
}

export type InputOutputConfig = {
  files?: string[]
  params?: string[]
}

export interface RunReturn {
  run_id: number
}

export interface RunStatus {
  [k: string]: StatusStringFine
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

export const STATUS = {
  PENDING: 'pending',
  UPCOMING: 'upcoming',
  COMPLETE: 'complete',
  ERROR: 'error',
  RUNNING: 'running'
} as const;

export type StatusString = typeof STATUS[keyof typeof STATUS];

export interface StepStatus {
  name: string;
  status: StatusString;
  statusFine: StatusStringFine;
  error?: string;
}

export const defaultStepStatus: StepStatus = {
  name: '',
  status: 'pending',
  statusFine: 'NOT STARTED',
};

export const defaultWorkspaceStatus = {
  steps: {} as {[k: string]: StepStatus},
  output_files: [] as string[],
  file_contents: {} as {[k: string]: any}, // key = filenames
  last_params: {} as {[k: string]: any},
  params: {} as {[k: string]: {[k: string]: Maybe<any>}}, // top key = 'param1/param2/...paramN' if many from 1; innner keys = param output's keys.
  ui_contents: {} as {[k: string]: {[k: string]: Maybe<any>}}, // top key = name defined in vulcan config (matches step name in dag); inner keys = file outputs' filenames
}

export type WorkspaceStatus = typeof defaultWorkspaceStatus;

// The elements stored as local cookie.
export type VulcanStorage = {
  workspace: Workspace,
  status: WorkspaceStatus
};

export const defaultVulcanStorage: VulcanStorage = {
  workspace: {...defaultWorkspace},
  status: {...defaultWorkspaceStatus},
};

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
