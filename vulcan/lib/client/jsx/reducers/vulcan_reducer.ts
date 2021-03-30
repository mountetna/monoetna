import {
    VulcanAction
} from '../actions/vulcan';
import {defaultSessionStatusResponse, SessionStatusResponse, Workflow, WorkflowsResponse} from "../api_types";


export type DownloadedData = any; // TODO: improve typing here.
export type DownloadedStepDataMap = {[k: string]: DownloadedData};

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
    calculating: false,
};

export type VulcanState = Readonly<(typeof defaultVulcanState)>;

// Copy over any downloaded data that is still referenced in the state.status download maps,
// freeing resources that are no longer necessary.
export function updateDownloads(state: VulcanState): VulcanState {
    let data: Readonly<DownloadedStepDataMap> = {};

    state.status[0].forEach(step => {
        for (let k in step.downloads) {
            const url = step.downloads[k];
            if (url) {
                data = {...data, [url]: state.data[url]};
            } else {
                console.warn('Unexpected empty url for download', k, step.downloads)
            }
        }
    })

    return {...state, data: data};
}

export default function VulcanReducer(state: VulcanState, action: VulcanAction): VulcanState {
    state = state || defaultVulcanState;

    switch (action.type) {
        case 'SET_WORKFLOWS':
            return {
                ...state,
                workflows: action.workflows,
            };
        case 'SET_WORKFLOW':
            return {
                ...state,
                workflow: action.workflow,
                session: {...state.session, workflow_name: action.workflow.name, key: ""},
            };
        case 'SET_STATUS':
            state = {...state, status: action.status};
            state = updateDownloads(state);
            return state;

        case 'SET_DOWNLOAD':
            return {
                ...state,
                data: {
                    ...state.data,
                    [action.url]: action.data
                }
            };

        case "RELEASE_DOWNLOAD":
            const data = { ...state.data };
            delete data[action.url];

            return {
                ...state,
                data
            };

        case 'SET_SESSION':
            return {
                ...state,
                session: { ...action.session, project_name: CONFIG.project_name },
            };

        case 'SET_INPUTS':
            return {
                ...state,
                inputs: action.inputs,
            }

        case 'COMMIT_INPUTS':
            // Ensure the reference is shared, making equality comparison simpler.
            const newInputs = state.inputs || {};
            return {
                ...state,
                inputs: newInputs,
                session: {
                    ...state.session,
                    inputs: newInputs,
                }
            }

        case 'SET_CALCULATING':
            return {
                ...state,
                calculating: action.calculating,
            };

        default:
            return state;
    }
}
