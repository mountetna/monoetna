import {VulcanAction} from '../actions/vulcan_actions';
import {
  Workflow,
  Workspace,
  AccountingReturn,
  defaultWorkspaceStatus,
  defaultWorkflow,
  StatusStringBroaden,
  WorkspaceMinimal
} from '../api_types';
import {
  allUIStepNames,
  inputUINames,
  paramValuesToRaw,
  uiContentsFromFiles,
  upcomingStepNames,
  vulcanConfigFromRaw,
} from '../selectors/workflow_selectors';
import {Maybe, withDefault} from '../selectors/maybe';
import { pick } from 'lodash';
import { stepOutputMapping } from '../components/workspace/ui_definitions/input_types';

const defaultWorkflows = [] as Workflow[];
const defaultWorkspaces = [] as WorkspaceMinimal[];
const defaultWorkspace = null as Workspace | null;
const defaultId = null as number | null;
const defaultStatus = defaultWorkspaceStatus;
const defaultValidationErrors: [string | null, string, string[]][] = [];

export const defaultVulcanState = {
  projectName: '',
  workflows: defaultWorkflows,
  workflow: defaultWorkflow,
  workspaces: defaultWorkspaces,
  workspace: defaultWorkspace,
  workspaceId: defaultId,
  configId: defaultId,
  runId: defaultId,
  isRunning: false,

  // rule/step statuses, ui&param contents per local, file&param contents per server 
  status: defaultStatus, // Only filled input/outputs in here

  // Trigger re-pull via useEffect inside vulcan_context
  update_workflows: true,
  update_files: false, // ui-related files needed for the workspace-manager

  // Trigger send to server
  pushSteps: [] as string[],

  // MAYBE: So that users can preview from their local storage before triggering the vulcan workspace to match it
  // Not used yet
  sessionWorkspaceSyncd: true,

  // A subset of all steps that have buffered changes.
  bufferedSteps: [] as (string | null)[],
  // If a buffered step was filled in while remaining 'pending'
  workQueueable: false,

  // Step marked for auto-passing by user
  autoPassSteps: [] as (string | null)[],

  // Request Run and trigger polling state
  triggerRun: [] as (string | null)[],
  pollingState: 0,
  
  validationErrors: defaultValidationErrors,
};

export type VulcanState = Readonly<typeof defaultVulcanState>;

function useAccounting(
  state: VulcanState,
  action: {type: 'USE_UI_ACCOUNTING'} & {
    accounting: AccountingReturn;
    submittingStep: string;
    removeSync: boolean;
  }
): VulcanState {
  /*
  v1: "filterStaleness" outputs matching hash of step inputs were used for determining staleness
  v2: snakemake determines staleness and the back-end returns:
        - 'config_id' = an id unique to the workspace setup
        - 'scheduled' = will run in next Run, and
        - 'downstream' = not ready to run yet, but is downstream in the dag of scheduled steps.
      Filters staleness, updates step statuses, and sets the new configId.
  */
  const {workspace, status} = state;
  if (!workspace) return state;

  let newStatus = {...status};
  // Use returned params as source of truth
  newStatus.last_params = action.accounting.params;
  // Clear ui_content knowledge and output_files knowledge. We will fill it back in next via 'update_files: true' below.
  newStatus.output_files = [];
  newStatus.ui_contents = uiContentsFromFiles(workspace);

  let newSteps = {...newStatus.steps};
  for (let step of action.accounting.scheduled.concat(action.accounting.downstream)) {
    // The submitting step is pushing a new value from the client up, thus
    // it should not have its input made stale.
    if (step === action.submittingStep || !Object.keys(newSteps).includes(step)) {
      console.log(`${step} not tracked, skipping`)
      continue;
    }

    // Update steps' statuses
    newSteps[step] = {
      name: step,
      statusFine: "NOT STARTED",
      status: action.accounting.scheduled.includes(step) ? 'upcoming' : 'pending'
    }
  }

  return {
    ...state,
    status: {...newStatus, steps: {...newSteps}},
    workQueueable: upcomingStepNames(workspace, newStatus).length > 0,
    configId: action.accounting.config_id,
    update_files: true,
    pushSteps: action.removeSync ? state.pushSteps.filter(s => s!=action.submittingStep) : state.pushSteps,
  };
}

export default function VulcanReducer(
  state: VulcanState,
  action: VulcanAction
): VulcanState {
  state = state || defaultVulcanState;

  switch (action.type) {
    // case 'SET_PROJECT':
    //   return {
    //     ...state,
    //     projectName: action.project
    //   };
    // case 'SET_WORKFLOWS':
    //   return {
    //     ...state,
    //     workflows: action.workflows
    //   };
    // case 'SET_WORKSPACES':
    //   return {
    //     ...state,
    //     workspaces: action.workspaces
    //   };
    case 'MODIFY_POLLING':
      return {
        ...state,
        pollingState: Math.max(state.pollingState + action.delta, 0)
      };
    case 'SET_WORKFLOW':
      const workflowProject = action.workflow.project_name;
      if (!["all", action.projectName].includes(workflowProject)) {
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
      // // ToDo: Ensure Project
      // const workspaceProject = action.workspace.project;
      // if (!workspaceProject) {
      //   // attempt to fill project
      //   const workflow = workflowByIdFromWorkflows(action.workspace.workflow_id, state.workflows);
      //   action.workspace.project = workflow?.project
      // } else if (workspaceProject!=action.projectName) {
      //   // Reject if somehow current project doesn't match the workspace project
      //   return state;
      // }

      return {
        ...state,
        workspace: {
          ...action.workspace,
          vulcan_config: vulcanConfigFromRaw(action.workspace.vulcan_config)
        }
      };
    case 'UPDATE_WORKFLOWS_AND_WORKSPACES':
      return {
        ...state,
        update_workflows: true
      }
    case 'SET_WORKFLOWS_AND_WORKSPACES':
      const updates = pick(action.updates, ['workflows', 'workspaces', 'workspace'])
      return {
        ...state,
        ...updates,
        update_workflows: false
      }
    case 'SET_FULL_WORKSPACE_STATUS':
      if (state.workspaceId !== action.workspace.workspace_id) {
        console.warn('Cannot update to show workspace of the wrong id.');
      }
      return {
        ...state,
        workspace: {...action.workspace},
        status: {...action.status},
        update_files: action.update_files,
        isRunning: action.isRunning,
        configId: action.workspace.last_config_id,
        runId: action.workspace.last_run_id,
      }
    case 'SET_CONFIG_ID':
      return {
        ...state,
        configId: action.configId
      };
    case 'SET_RUN_ID':
      return {
        ...state,
        runId: action.runId
      };
    case 'SET_LAST_CONFIG':
      return {
        ...state,
        status: {
          ...state.status,
          last_params: action.lastConfig
        }
      };
    case 'USE_UI_ACCOUNTING':
      // Arrive here from sending ui-steup / config to the back-end
      // When a submitting step is given, filter stale inputs that result from submitting a change
      // to that step's outputs in session.inputs.
      // ToDo once cache'ing: Also assess if 'stale' inputs can be filled in with versions matching sent setup.
      return useAccounting(state, action);
    case 'SET_STATUS_FROM_STATUSES':
      if (!action.statusReturns) {
        console.log("empty data")
      }
      // Arrive here from polling return
      const newStepStatus = {};
      const newCompletions = [] as string[];
      for (let [stepName, statusFine] of Object.entries(action.statusReturns)) {
        const statusBroad = StatusStringBroaden(statusFine);
        newStepStatus[stepName] = {
          name: stepName,
          status: statusBroad,
          statusFine: statusFine
        }
        if (statusBroad=='complete' && state.status.steps[stepName].status!='complete') {
          newCompletions.push(stepName);
        }
      };
      const newStatus = {
        ...state.status,
        steps: {...state.status.steps, ...newStepStatus}
      };
      return {
        ...state,
        isRunning: action.isRunning,
        status: {...newStatus},
        workQueueable: upcomingStepNames(state.workspace as Workspace, newStatus).length > 0,
        update_files: newCompletions.length > 0
      };

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

    case 'SET_STATE_FROM_STORAGE':
      const currentWorkspace = state.workspace;
      if (
        !currentWorkspace ||
        action.storage.workspace.workflow_id !== currentWorkspace.workflow_id ||
        action.storage.workspace.workspace_id !== currentWorkspace.workspace_id //||
        // action.storage.workspace.project !== currentWorkspace.project
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

      if (!state.workspace || !state.workspace.vulcan_config || !(action.stepName in state.workspace.vulcan_config)) {
        return state;
      }

      // Map elements from ui-intenal keys to external keys (param and filenames)
      const newValues = {}
      const mapping = stepOutputMapping(state.workspace.vulcan_config[action.stepName])
      Object.entries(mapping).forEach(
        ([internal, external]) => {
          newValues[external] = action.values[internal]
        }
      )

      // Set in status, and mark for sending to the server
      const statusElement = allUIStepNames(state).includes(action.stepName) ? 'ui_contents' : 'params';
      return {
        ...state,
        pushSteps: state.pushSteps.concat(action.stepName),
        status: {
          ...state.status,
          [statusElement]: {
            ...state.status[statusElement],
            [action.stepName]: newValues
          },
        }
      };

    case 'REMOVE_SYNC':
      return {
        ...state,
        pushSteps: state.pushSteps.filter(s => s!=action.stepName)
      };

    case 'UPDATE_FILES':
      return {
        ...state,
        status: {
          ...state.status,
          ...action.statusUpdates,
          ui_contents: uiContentsFromFiles(state.workspace as Workspace, action.statusUpdates.file_contents)
        },
        update_files: false
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
      if (action.step == null) {
        ready = false;
      } else {
        const stepNum = state.status[0].findIndex(
          (element) => element['name'] == stepName
        );
        ready = state.status[0][stepNum]['status'] == 'pending';
      }
      return {
        ...state,
        workQueueable: ready
      };

    case 'CLEAR_RUNNING':
      return {
        ...state,
        isRunning: false
      }

    case 'CLEAR_CHANGES_READY':
      return {
        ...state,
        workQueueable: false
      };

    default:
      return state;
  }
}
