import {
  SET_DATA,
  SET_WORKFLOW,
  SET_WORKFLOWS,
  SET_STATUS,
  SET_PATH,
  SET_STEP
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
      // Status returned by the server as an Array of Arrays,
      //   with possible looping paths as the nested Arrays.
      // Assume stable order from the server.
      let statusWorkflow = {...state.workflow};

      statusWorkflow.steps.forEach((path, pathIndex) => {
        path.forEach((step, stepIndex) => {
          statusWorkflow.steps[pathIndex][stepIndex] = {
            ...step,
            ...action.status[pathIndex][stepIndex]
          };
        });
      });
      return {
        ...state,
        workflow: statusWorkflow
      };
    case SET_DATA:
      // Inject the data payload to state based on the URL.
      // Workflow steps are an Array of Arrays, so need to
      //   loop through them to find matching data URLs.
      // Assume stable order from the server.
      let dataWorkflow = {...state.workflow};

      dataWorkflow.steps.forEach((path, pathIndex) => {
        const stepIndex = path.findIndex((step) => {
          return step.data_url === action.url;
        });

        dataWorkflow.steps[pathIndex][stepIndex] = {
          ...dataWorkflow.steps[pathIndex][stepIndex],
          data: action.data
        };
      });

      return {
        ...state,
        workflow: dataWorkflow
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
    default:
      return state;
  }
}
