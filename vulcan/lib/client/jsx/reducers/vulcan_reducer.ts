import {
  VulcanAction
} from '../actions/vulcan_actions';
import {
  defaultSessionStatusResponse,
  SessionStatusResponse, StepStatus,
  Workflow,
  WorkflowsResponse
} from "../api_types";
import {allExpectedOutputSources, filterEmptyValues, statusOfStep, stepOfStatus} from "../selectors/workflow_selectors";
import {mapSome, Maybe} from "../selectors/maybe";

export type DownloadedData = any; // TODO: improve typing here.
export type DownloadedStepDataMap = { [k: string]: DownloadedData };

const defaultWorkflows: WorkflowsResponse['workflows'] = [];
const defaultWorkflow: Workflow | null = null;
const defaultStatus: SessionStatusResponse['status'] = [[]];
const defaultData: DownloadedStepDataMap = {};
const defaultInputs: SessionStatusResponse['session']['inputs'] = {};
export const defaultSession: SessionStatusResponse['session'] = {
  project_name: '',
  workflow_name: '',
  key: '',
  inputs: {},
}
const defaultValidationErrors: [string | null, string, string[]][] = [];

export const defaultVulcanState = {
  workflows: defaultWorkflows,
  workflow: defaultWorkflow as Workflow | null,
  status: defaultStatus,
  data: defaultData,

  // a subset of all inputs that might have changes to them.
  bufferedInputValues: {} as {[k: string]: Maybe<any>},
  // a subset of all inputs that have been explicitly edited by the user,
  // and the hash reported by status for that edit (if it exists).
  // used to clear inputs /outputs that are stale, such as when an input
  // change has downstream effects to the decision of an input.
  stepChangedAt: {} as {[k: string]: Maybe<string>},

  session: defaultSession,
  outputs: defaultSessionStatusResponse.outputs,
  validationErrors: defaultValidationErrors,
  pollingState: 0,
  editedSource: null as null | string,
};

export type VulcanState = Readonly<(typeof defaultVulcanState)>;

export default function VulcanReducer(state: VulcanState, action: VulcanAction): VulcanState {
  state = state || defaultVulcanState;

  switch (action.type) {
    case 'SET_WORKFLOWS':
      return {
        ...state,
        workflows: action.workflows,
      };
    case 'MODIFY_POLLING':
      return {
        ...state,
        pollingState: Math.max(state.pollingState + action.delta, 0)
      };
    case 'SET_WORKFLOW':
      const workflowProjects = action.workflow.projects;
      if (!workflowProjects || !workflowProjects.includes(action.projectName)) return state;

      return {
        ...state,
        workflow: action.workflow,
        data: defaultData,
        session: {...defaultSession, workflow_name: action.workflow.name, project_name: action.projectName},
      };
    case 'SET_STATUS':
      const newInputs = filterStaleInputs(state, action);
      const newData = filterStaleData(state, action);

      return {
        ...state,
        stepChangedAt: {},
        status: action.status,
        session: {...state.session, inputs: newInputs},
        data: newData,
      };

    case 'SET_BUFFERED_INPUT':
      return {...state, bufferedInputValues: action.inputs};

    case 'CLEAR_BUFFERED_INPUT':
      return {...state, bufferedInputValues: {}};

    case 'SET_DOWNLOAD':
      return {
        ...state,
        data: {
          ...state.data,
          [action.url]: action.data
        }
      };

    case 'SET_SESSION':
      // Ignore sessions not for this workflow, project combination. safety guard.
      // The project name and workflow names are locked in via the set workflow action.
      const sessionWorkflow = state.workflow;
      if (!sessionWorkflow || action.session.workflow_name !== sessionWorkflow.name || action.session.project_name !== state.session.project_name) {
        console.error('Cannot set session, project name / workflow name does not match.')
        return state;
      }

      return {
        ...state,
        session: {...state.session, ...action.session},
      };

    case 'REMOVE_DOWNLOADS':
      return {
        ...state,
        status: [
          state.status[0].map(status => action.stepNames.includes(status.name) ? {...status, downloads: null} : status),
        ],
      };

    case "REMOVE_INPUTS":
      if (state.pollingState) {
        console.error('cannot change inputs while polling...')
        return state;
      }

      const removedInputs = {...state.session.inputs};
      action.inputs.forEach(source => delete removedInputs[source]);

      return {
        ...state,
        session: {...state.session, inputs: filterEmptyValues(removedInputs)},
      }


    case 'SET_INPUTS':
      if (state.pollingState) {
        console.error('cannot change inputs while polling...')
        return state;
      }

      return {
        ...state,
        session: {...state.session, inputs: filterEmptyValues(action.inputs)},
      }

    case 'ADD_VALIDATION_ERRORS':
      return {
        ...state,
        validationErrors: [
          ...state.validationErrors,
          [action.stepName, action.inputLabel, action.errors]
        ]
      }

    case 'REMOVE_VALIDATION_ERRORS':
      return {
        ...state,
        validationErrors: state.validationErrors.filter(([_1, _2, e]) => e !== action.errors)
      }

    default:
      return state;
  }
}

function filterStaleInputs(state: VulcanState, action: { type: "SET_STATUS" } & { status: [StepStatus[]] }) {
  let newInputs = {...state.session.inputs};
  for (let stepName in state.stepChangedAt) {
    mapSome(state.stepChangedAt[stepName], expectedHash => {
      if (!state.workflow) return;
      const step = stepOfStatus(stepName, state.workflow);
      if (!step) return;
      const newStatus = statusOfStep(stepName, action.status);
      if (newStatus && newStatus.hash === expectedHash) return;

      allExpectedOutputSources(step).forEach(outputSource => {
        delete newInputs[outputSource];
      })
    })
  }
  return newInputs;
}

function filterStaleData(state: VulcanState, action: { type: "SET_STATUS" } & { status: [StepStatus[]] }) {
  let newData = {...defaultData};
  for (let oldStatus of state.status[0]) {
    const newStatus = statusOfStep(oldStatus.name, action.status);
    if (newStatus && newStatus.hash === oldStatus.hash && newStatus.downloads) {
      Object.values(newStatus.downloads).forEach(url => {
        if (url in state.data) {
          newData[url] = state.data[url];
        }
      })
    }
  }

  return newData;
}
