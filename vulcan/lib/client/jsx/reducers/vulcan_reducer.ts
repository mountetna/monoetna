import {
  VulcanAction
} from '../actions/vulcan';
import {defaultSessionStatusResponse, SessionStatusResponse, Workflow, WorkflowsResponse} from "../api_types";
import {unsetDependentInputs, filterEmptyValues} from "../selectors/workflow_selectors";


export type DownloadedData = any; // TODO: improve typing here.
export type DownloadedStepDataMap = { [k: string]: DownloadedData };

const defaultWorkflows: WorkflowsResponse['workflows'] = [];
const defaultWorkflow: Workflow | null = null;
const defaultStatus: SessionStatusResponse['status'] = [[]];
const defaultData: DownloadedStepDataMap = {};
const defaultInputs: SessionStatusResponse['session']['inputs'] = {};
const defaultSession: SessionStatusResponse['session'] = {
  project_name: '',
  workflow_name: '',
  key: '',
  inputs: {},
}

export const defaultVulcanState = {
  workflows: defaultWorkflows,
  workflow: defaultWorkflow as Workflow | null,
  status: defaultStatus,
  data: defaultData,
  inputs: defaultInputs,
  session: defaultSession,
  outputs: defaultSessionStatusResponse.outputs,
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
    case 'SET_WORKFLOW':
      const workflowProjects = action.workflow.projects;
      if (!workflowProjects || workflowProjects.length < 1) return state;

      return {
        ...state,
        workflow: action.workflow,
        session: {...defaultSession, workflow_name: action.workflow.name, project_name: workflowProjects[0]},
      };
    case 'SET_STATUS':
      return {...state, status: action.status};

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
      if (!sessionWorkflow) return state;
      if (action.session.workflow_name !== sessionWorkflow.name) return state;
      const projects = sessionWorkflow.projects;
      if (!projects || projects.length < 1) return state;
      if (action.session.project_name !== projects[0]) return state;

      return {
        ...state,
        session: action.session,
        inputs: {...state.inputs, ...action.session.inputs},
      };

    case 'REMOVE_DOWNLOADS':
      return {
        ...state,
        status: [
          state.status[0].map(status => action.stepNames.includes(status.name) ? {...status, downloads: null} : status),
        ],
      };

    case "REMOVE_INPUTS":
      const removedInputs = {...state.inputs};
      action.inputs.forEach(source => delete removedInputs[source]);

      return {
        ...state,
        inputs: removedInputs,
        session: {...state.session, inputs: filterEmptyValues(removedInputs)},
      }


    case 'SET_INPUTS':
      return {
        ...state,
        inputs: action.inputs,
        session: {...state.session, inputs: filterEmptyValues(action.inputs)},
      }

    case 'PATCH_INPUTS':
      const updatedInputs = unsetDependentInputs(action.inputs, {...state.inputs, ...action.inputs}, state.workflow);

      return {
        ...state,
        inputs: updatedInputs,
        session: {...state.session, inputs: filterEmptyValues(updatedInputs)},
      }

    default:
      return state;
  }
}
