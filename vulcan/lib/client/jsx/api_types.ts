import { Maybe } from "./selectors/maybe";

// CWL Step RUN Sentinels
export const RUN = {
  UI_QUERY: 'ui-queries/',
  UI_OUTPUT: 'ui-outputs/'
};

export const defaultWorkflowsResponse = {} as {workflows: Workflow[]};

export type WorkflowsResponse = typeof defaultWorkflowsResponse;
// export type ProjectWorkflowsResponse = typeof defaultWorkflowsResponse;
export type WorkflowCreateResponse = {'workflow_id': number, 'workflow_name': string};

// export interface StepInput {
//   id: string;
//   source: string;
//   doc?: string | null;
//   label?: string | null;
// }

export interface WorkspaceStep {
  name: string;
  input: InputConfig;
  output: OutputConfig;
  label?: string;
  ui_component?: string;
  doc?: string;
  // for stepInputs where a default may come in from the vulcan_config
  default?: any;
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

export type WorkflowServer = {
  id: number,
  project_name: string,
  name: string,
  repo_remote_url: string,
  created_at: number,
  updated_at: number,
}

export interface Workflow {
  id: number | null; // null to indicate stubbage
  name: string;
  project_name: string;
  repo_remote_url: string;
  created_at: string;
  updated_at: string;
  vignette?: string;
  image?: string;
}

export const defaultWorkflow: Workflow = {
  id: null,
  name: '',
  project_name: '',
  repo_remote_url: '',
  created_at: '0',
  updated_at: '0',
};

export type CreateWorkspaceResponse = {
  workspace_id: number,
  workflow_id: number,
  vulcan_config: VulcanConfigRaw,
  dag: string[]
}


// Per getWorkspace*s*
export interface WorkspaceMinimalMinusInconsistent {
  workspace_id: number | null;
  name: string;
  workflow_name: string;
  workflow_id: number | null;
  user_email: string;
  tags: string[] | null;
  git_ref: string;
  git_sha: string;
  dag: string[];
  created_at: string;
  updated_at: string;
  workspace_path: string;
  // NEEDED
  thumbnails: any[]
}
// ToDo: Remove these if the final version is stably not edited!
// Not currently in the minimal: none
export interface WorkspaceMinimalRaw extends WorkspaceMinimalMinusInconsistent {}
export interface WorkspaceMinimal extends WorkspaceMinimalMinusInconsistent {}
export type WorkspacesResponseRaw = {workspaces: WorkspaceMinimalRaw[]}
export type WorkspacesResponse = {workspaces: WorkspaceMinimal[]}

// Per getWorkspace
interface WorkspaceMinusInconsistent {
  workspace_id: number | null;
  workflow_id: number | null;
  workflow_name: string;
  name: string;
  user_email: string;
  workspace_path: string;
  dag: string[];
  git_ref: string;
  git_sha: string;
  created_at: string;
  updated_at: string;
  last_config: {[k: string]: any} | null;
  last_config_id: number | null;
  last_job_status: {[k: string]: StatusStringFine} | null;
  last_run_id: number | null;
  tags: string[] | null;
}
export interface WorkspaceRaw extends WorkspaceMinusInconsistent {
  vulcan_config: VulcanConfigElement[];
  // target_mapping: {[k: string]: WorkspaceStep};
}
export interface Workspace extends WorkspaceMinusInconsistent {
  vulcan_config: VulcanConfig;
  // project?: string;
  vignette?: string;
  thumbnails?: string[];
}

export const defaultWorkspace: Workspace = {
  workspace_id: null,
  workflow_id: null,
  workflow_name: '',
  name: '',
  user_email: '',
  workspace_path: '',
  tags: [],
  dag: [],
  git_ref: 'main',
  git_sha: '',
  vulcan_config: {},
  created_at: '',
  updated_at: '',
  last_config: {},
  last_config_id: null,
  last_job_status: {},
  last_run_id: null,
};

export interface FileContentResponse {
  filename: string;
  content: string;
  encoding: string;
}

export type MultiFileContentResponse = {files: FileContentResponse[]}

export interface MultiFileContent {
  [filename: string]: any;
}

export type FlatParams = {
  [k: string]: any
};

export interface AccountingReturn {
  config_id: number,
  params: {[k: string]: any},
  scheduled: string[],
  downstream: string[],
}

export interface isRunningReturn {
  running: boolean
}

// export type VulcanConfig = VulcanConfigElement[]
export type VulcanConfig = {[k: string]: VulcanConfigElement}
export type VulcanConfigRaw = VulcanConfigElement[];

export type VulcanConfigElement = {
  name: string;
  display: string;
  ui_component: string;
  label?: string;
  default?: any;
  doc?: string;
  input?: InputConfig;
  output?: OutputConfig;
  await_files?: string[];
}

export type InputConfig = {
  files?: string[] | {[k: string]: string}
  params?: string[] | {[k: string]: string}
  preset?: {[k: string]: any}
}

export type OutputConfig = {
  files?: string[] | {[k: string]: string}
  params?: string[] | {[k: string]: string}
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

export type LatencyReturn = {
  latency:  string // of form "#{median_latency}ms"
};

export type ClusterStatusReturn = {
  connection_success: boolean;
  expected_down: boolean;
  message: string; // of form "#{median_latency}ms"
};