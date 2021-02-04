import {
  SET_DATA,
  SET_WORKFLOW,
  SET_WORKFLOWS,
  SET_STATUS
} from '../actions/archimedes';

export default function ArchimedesReducer(state, action) {
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
      // Inputs probably returned by the server as Hash.
      // Need to re-map them to the workflow inputs, which
      //   should be an Array.
      let statusWorkflow = {...state.workflow};

      statusWorkflow.steps.forEach((step, index) => {
        statusWorkflow.steps[index] = {
          ...step,
          ...action.status[step.name]
        };
      });
      return {
        ...state,
        workflow: statusWorkflow
      };
    case SET_DATA:
      // Inject the data to state based on the URL.
      let dataWorkflow = {...state.workflow};
      const stepIndex = dataWorkflow.steps.findIndex((step) => {
        return step.data_url === action.url;
      });
      dataWorkflow.steps[stepIndex] = {
        ...dataWorkflow.steps[stepIndex],
        data: action.data
      };

      return {
        ...state,
        workflow: dataWorkflow
      };
    default:
      return state;
  }
}
