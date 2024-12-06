import {VulcanAction} from '../actions/vulcan_actions';
import {
  WorkflowsResponse,
  Workspace,
  AccountingReturn,
  defaultWorkspaceStatus,
  WorkspaceStatus,
  Workspaces,
  defaultWorkflow
} from '../api_types';
import {
  allUISteps,
  configUISteps,
} from '../selectors/workflow_selectors';
import {mapSome, Maybe, some, withDefault} from '../selectors/maybe';

export type DownloadedData = any; // TODO: improve typing here.
export type DownloadedStepDataMap = {[k: string]: DownloadedData};

const defaultWorkflows = [] as WorkflowsResponse;
const defaultWorkspaces = [] as Workspaces;
const defaultWorkspace = null as Workspace | null;
const defaultId = null as number | null;
const defaultStatus = defaultWorkspaceStatus;
const defaultFiles = [] as string[];
const defaultData: DownloadedStepDataMap = {};
// const defaultInputs: SessionStatusResponse['session']['inputs'] = {};
// export const defaultSession: SessionStatusResponse['session'] = {
//   project_name: '',
//   workflow_id: null,
//   key: '',
//   inputs: {},
//   reference_figure_id: null
// };
// const defaultFigure: VulcanFigure = {
//   figure_id: null,
//   inputs: {},
//   id: null
// };
const defaultValidationErrors: [string | null, string, string[]][] = [];

export const defaultVulcanState = {
  projectName: '',
  workflows: defaultWorkflows,
  workflow: defaultWorkflow,
  workspaces: defaultWorkspaces,
  workspaceId: defaultId,
  workspace: defaultWorkspace, // All possible input/outputs defined in here
  configId: defaultId,
  runId: defaultId,
  status: defaultStatus, // Only filled input/outputs in here
  data: defaultData,

  // MAYBE: So that users can preview from their local storage before triggering the vulcan workspace to match it
  sessionWorkspaceSyncd: true,

  // A subset of all steps that have buffered changes.
  bufferedSteps: [] as (string | null)[],
  // If a buffered step was filled in while remaining 'pending'
  committedStepPending: false,

  // Step marked for auto-passing by user
  autoPassSteps: [] as (string | null)[],
  triggerRun: [] as (string | null)[],
  
  validationErrors: defaultValidationErrors,
  pollingState: 0
};

export type VulcanState = Readonly<typeof defaultVulcanState>;

function filterStaleness(
  state: VulcanState,
  action: {type: 'SET_STATUS'} & {
    accounting: AccountingReturn;
    submittingStep: Maybe<string>;
  }
) {
  /*
  v1: outputs matching hash of step inputs were used for determining staleness
  v2: snakemake determines staleness and the back-end returns 'scheduled' = will run in next Run, and
      'downstream' = not ready to run yet, but is downstream in the dag of scheduled steps
  */
  const {workspace, status, data} = state;
  if (!workspace) return state;
  
  let newStatus = {...status};
  let newData = {...data}

  const staleSteps = action.accounting.downstream.concat(action.accounting.scheduled)
  const staleUISteps = staleSteps.filter(name => allUISteps(state).includes(name))

  return withDefault(
    mapSome(action.submittingStep, (submittingStep) => {

      for (let [step, stepStatus] of Object.entries(newStatus.steps)) {
        // The submitting step is pushing a new value from the client up, thus
        // it should not have its input made stale.
        if (step === submittingStep || !staleSteps.includes(step)) continue;

        // Clear ui_content knowledge, output_files knowledge, and downloaded output file content
        // (May still exist in workspace, but we will assume it's stale.)
        if (stepStatus.outputs?.files) {
          for (let output in stepStatus.outputs.files) delete newData[output];
          newStatus.output_files.filter(name => !stepStatus.outputs?.files?.includes(name));
          if (staleUISteps.includes(step)) {
            delete newStatus.ui_contents[step];
          }
        }

        // Update steps' statuses
        delete newStatus.steps[step].error ;
        delete newStatus.steps[step].outputs;
        newStatus.steps[step].statusFine = "NOT STARTED";
        newStatus.steps[step].status = action.accounting.scheduled.includes(step) ? "upcoming" : "pending";
      }
      
      return {...state, status: newStatus};
    }),
    state
  );
}

export default function VulcanReducer(
  state: VulcanState,
  action: VulcanAction
): VulcanState {
  state = state || defaultVulcanState;

  switch (action.type) {
    case 'SET_PROJECT':
      return {
        ...state,
        projectName: action.project
      };
    case 'SET_WORKFLOWS':
      return {
        ...state,
        workflows: action.workflows
      };
    case 'SET_WORKSPACES':
      return {
        ...state,
        workspaces: action.workspaces
      };
    case 'MODIFY_POLLING':
      return {
        ...state,
        pollingState: Math.max(state.pollingState + action.delta, 0)
      };
    case 'SET_WORKFLOW':
      const workflowProjects = action.workflow.projects;
      if (
        !workflowProjects.includes("all") && // CHECK ME!!!!!
          !workflowProjects.includes(action.projectName)
      ) {
        return state;
      }

      return {
        ...state,
        workflow: action.workflow
      };
    case 'SET_WORKSPACE_ID':
      return {
        ...state,
        workspaceId: action.workspaceId
      };
    case 'SET_WORKSPACE':
      const workspaceProject = action.workspace.project;
      if (workspaceProject!=action.projectName) {
        return state;
      }

      return {
        ...state,
        workspace: action.workspace
      };
    case 'SET_CONFIG_ID':
      return {
        ...state,
        workspaceId: action.configId
      };
    case 'SET_RUN_ID':
      return {
        ...state,
        workspaceId: action.runId
      };
    case 'SET_STATUS':
      // When a submitting step is given, filter stale inputs that result from submitting a change
      // to that step's outputs in session.inputs.
      state = filterStaleness(state, action);
      return {...state};

    case 'SET_BUFFERED_INPUT':
      if (state.bufferedSteps.includes(action.step)) {
        return state;
      }
      return {...state, bufferedSteps: [...state.bufferedSteps, action.step]};

    case 'CLEAR_BUFFERED_INPUT':
      const idx = state.bufferedSteps.indexOf(action.step);
      if (idx === -1) return state;
      const bufferedSteps = [...state.bufferedSteps];
      bufferedSteps.splice(idx, 1);
      return {...state, bufferedSteps};

    case 'SET_AUTO_PASS_STEP':
      if (state.autoPassSteps.includes(action.step)) {
        return state;
      }
      return {...state, autoPassSteps: [...state.autoPassSteps, action.step]};

    case 'CLEAR_AUTO_PASS_STEP':
      const index = state.autoPassSteps.indexOf(action.step);
      if (index === -1) return state;
      const autoPassSteps = [...state.autoPassSteps];
      autoPassSteps.splice(index, 1);
      return {...state, autoPassSteps};

    case 'SET_RUN_TRIGGER':
      if (state.triggerRun.includes(action.step)) {
        return state;
      }
      return {...state, triggerRun: [...state.triggerRun, action.step]};

    case 'CLEAR_RUN_TRIGGERS':
      const triggerRun = [...state.triggerRun].filter(
        (val) => !action.steps.includes(val)
      );
      return {...state, triggerRun};

    case 'SET_DOWNLOAD':
      return {
        ...state,
        data: {
          ...state.data,
          [action.fileName]: action.fileData
        }
      };

    case 'SET_STATE_FROM_STORAGE':
      const currentWorkspace = state.workspace;
      if (
        !currentWorkspace ||
        action.storage.workspace.workflow_id !== currentWorkspace.workflow_id ||
        action.storage.workspace.workspace_id !== currentWorkspace.workspace_id ||
        action.storage.workspace.project !== currentWorkspace.project
      ) {
        console.warn(
          'Cannot update session from cookie: project name, workflow id, or workspace id do not match.',
          currentWorkspace,
          action
      );
        return state;
      }
      return {
        ...state,
        workspace: {...action.storage.workspace},
        status: {...action.storage.status},
        sessionWorkspaceSyncd: false,
      };

    case 'SET_WORKSPACE_STATE_SYNCD':
      return {
        ...state,
        sessionWorkspaceSyncd: true,
      };

    // case 'REMOVE_DOWNLOADS':
    //   return {
    //     ...state,
    //     status: [
    //       state.status[0].map((status) =>
    //         action.stepNames.includes(status.name)
    //           ? {...status, downloads: null}
    //           : status
    //       )
    //     ]
    //   };

    case 'SET_WORKSPACE_FILES':
      return {
        ...state,
        status: {...state.status, output_files: action.fileNames}
      };

    case 'SET_UI_VALUES':
      if (state.pollingState) {
        console.error('cannot change inputs while polling...');
        return state;
      }

      let config_values: WorkspaceStatus['ui_contents'] = {};
      let ui_values: WorkspaceStatus['ui_contents'] = {};
      for (let step in Object.keys(action.values)) {
        if (configUISteps(state).includes(step)) {
          config_values[step] = action.values[step];
        } else {
          ui_values[step] = action.values[step];
        }
      }

      return {
        ...state,
        status: {...state.status, config_contents: config_values, ui_contents: ui_values}
      };

    case 'ADD_VALIDATION_ERRORS':
      return {
        ...state,
        validationErrors: [
          ...state.validationErrors,
          [action.stepName, action.inputLabel, action.errors]
        ]
      };

    case 'REMOVE_VALIDATION_ERRORS':
      return {
        ...state,
        validationErrors: state.validationErrors.filter(
          ([_1, _2, e]) => e !== action.errors
        )
      };

    case 'CHECK_CHANGES_READY':
      const stepName = withDefault(action.step, null);
      let ready: boolean;
      if (stepName == null) {
        ready = false;
      } else {
        const stepNum = state.status[0].findIndex(
          (element) => element['name'] == stepName
        );
        ready = state.status[0][stepNum]['status'] == 'pending';
      }
      return {
        ...state,
        committedStepPending: ready
      };

    case 'CLEAR_CHANGES_READY':
      return {
        ...state,
        committedStepPending: false
      };

    default:
      return state;
  }
}
