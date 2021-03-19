import {
    SET_DATA,
    SET_WORKFLOW,
    SET_WORKFLOWS,
    SET_STATUS,
    SET_PATH,
    SET_STEP,
    SET_SESSION,
    SET_INPUTS,
    COMMIT_INPUTS, VulcanAction
} from '../actions/vulcan';
import {Workflow, WorkflowsResponse} from "../api/types";

export const defaultVulcanState = {
    workflows: {} as WorkflowsResponse['workflows'],
    workflow: null as (Workflow | null),

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
            return {
                ...state,
                workflow: action.workflow
            };
        case SET_STATUS:
            // Update each entry in the status Arrays.
            // Workflow status are an Array of Arrays, so need to
            //   loop through them.
            // Assume stable order from the server.
            let currentStatus = [...(state.status || [])];
            action.status.forEach((path, pathIndex) => {
                if (!currentStatus[pathIndex]) currentStatus[pathIndex] = [];
                path.forEach((step, stepIndex) => {
                    if (!currentStatus[pathIndex][stepIndex])
                        currentStatus[pathIndex][stepIndex] = {};

                    currentStatus[pathIndex][stepIndex] = {
                        ...currentStatus[pathIndex][stepIndex],
                        ...step
                    };
                });
            });
            return {
                ...state,
                status: currentStatus
            };
        // case SET_DATA:
        //     // Inject the data payload to status based on the URL.
        //     // Workflow status are an Array of Arrays, so need to
        //     //   loop through them to find matching data URLs.
        //     // Assume stable order from the server.
        //     state = {...state};
        //     const dataCache = state.dataCache = Object.assign({}, state.dataCache);
        //     dataCache[action.url] = action.data;
        //     const dataStatus = state.status = [...state.status];
        //
        //     dataStatus.forEach((path, pathIndex) => {
        //         path.forEach((step, stepIndex) => {
        //             if (!step.downloads) return;
        //
        //             Object.keys(step.downloads).forEach((downloadKey) => {
        //                 if (action.url === step.downloads[downloadKey]) {
        //                     if (!dataStatus[pathIndex][stepIndex].data) {
        //                         dataStatus[pathIndex][stepIndex].data = {};
        //                     }
        //
        //                     dataStatus[pathIndex][stepIndex].data[downloadKey] = action.data;
        //                 }
        //             });
        //         });
        //     });
        //
        //     return state;
        // case SET_PATH:
        //     return {
        //         ...state,
        //         pathIndex: action.pathIndex
        //     };
        // case SET_STEP:
        //     return {
        //         ...state,
        //         stepIndex: action.stepIndex
        //     };
        // case SET_SESSION:
        //     return {
        //         ...state,
        //         session: action.session
        //     };
        // case SET_INPUTS:
        //     return {
        //         ...state,
        //         inputs: {
        //             ...(state.inputs || {}),
        //             ...action.inputs
        //         },
        //     };
        // case COMMIT_INPUTS:
        //     // Ensure the reference is shared, making equality comparison simpler.
        //     const newInputs = state.inputs || {};
        //
        //     return {
        //         ...state,
        //         inputs: newInputs,
        //         session: {
        //             inputs: newInputs,
        //         }
        //     };
        default:
            return state;
    }
}
