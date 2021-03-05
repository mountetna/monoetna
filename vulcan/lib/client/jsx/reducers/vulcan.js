import {
  SET_DATA,
  SET_WORKFLOW,
  SET_WORKFLOWS,
  SET_STATUS,
  SET_PATH,
  SET_STEP,
  SET_SESSION,
  SET_INPUTS,
  SET_CALCULATING
} from '../actions/vulcan';

export default function VulcanReducer(state, action) {
  state = state ? state : {};

  switch (action.type) {
    case SET_WORKFLOWS:
      return {
        ...state,
        workflows: action.workflows
      };
    case SET_WORKFLOW:
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
    case SET_DATA:
      // Inject the data payload to status based on the URL.
      // Workflow status are an Array of Arrays, so need to
      //   loop through them to find matching data URLs.
      // Assume stable order from the server.
      let dataStatus = [...state.status];
      dataStatus.forEach((path, pathIndex) => {
        path.forEach((step, stepIndex) => {
          if (!step.downloads) return;

          Object.keys(step.downloads).forEach((downloadKey) => {
            if (action.url === step.downloads[downloadKey]) {
              if (!dataStatus[pathIndex][stepIndex].data) {
                dataStatus[pathIndex][stepIndex].data = {};
              }

              dataStatus[pathIndex][stepIndex].data[downloadKey] = action.data;
            }
          });
        });
      });

      return {
        ...state,
        status: dataStatus
      };
    case SET_PATH:
      return {
        ...state,
        pathIndex: action.pathIndex
      };
    case SET_STEP:
      return {
        ...state,
        stepIndex: action.stepIndex
      };
    case SET_SESSION:
      return {
        ...state,
        session: action.session
      };
    case SET_INPUTS:
      console.log('state', state);
      console.log('action', action);
      return {
        ...state,
        session: {
          ...(state.session || {}),
          inputs: {
            ...(state.session ? state.session.inputs : {}),
            ...action.inputs
          }
        }
      };
    case SET_CALCULATING:
      // This should go away once we have async jobs and
      //   a different way to check status.
      return {
        ...state,
        calculating: action.calculating
      };
    default:
      return state;
  }
}
