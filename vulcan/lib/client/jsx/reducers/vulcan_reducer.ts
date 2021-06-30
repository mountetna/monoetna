import {
  VulcanAction
} from '../actions/vulcan_actions';
import {defaultSessionStatusResponse, SessionStatusResponse, Workflow, WorkflowsResponse} from "../api_types";
import {unsetDependentInputs, filterEmptyValues} from "../selectors/workflow_selectors";
import {Maybe, some} from "../selectors/maybe";

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
const defaultValidationErrors: [string, string[]][] = [];

export const defaultVulcanState = {
  workflows: defaultWorkflows,
  workflow: defaultWorkflow as Workflow | null,
  status: defaultStatus,
  data: defaultData,

  // When null, no step is being edited.
  // when some => null, the primary inputs are being edited.
  curEditingInputStep: null as Maybe<string | null>,

  // Maps output names to possibly values.
  bufferedInputValues: {} as {[k: string]: Maybe<any>},

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
        session: {...defaultSession, workflow_name: action.workflow.name, project_name: action.projectName},
      };
    case 'SET_STATUS':
      return {...state, status: action.status};

    case 'SET_BUFFERED_INPUT':
      return {...state, bufferedInputValues: action.value, curEditingInputStep: some(action.step)};

    case 'CLEAR_BUFFERED_INPUT':
      return {...state, bufferedInputValues: {}, curEditingInputStep: null};

    case 'SET_DOWNLOAD':
      return {
        ...state,
        data: {
          ...state.data,
          [action.url]: action.data
        }
      };

    case "RELEASE_DOWNLOAD":
      const data = {...state.data};
      delete data[action.url];

      return {
        ...state,
        data
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
          [action.inputLabel, action.errors]
        ]
      }

    case 'REMOVE_VALIDATION_ERRORS':
      return {
        ...state,
        validationErrors: state.validationErrors.filter(([_, e]) => e !== action.errors)
      }

    default:
      return state;
  }
}
