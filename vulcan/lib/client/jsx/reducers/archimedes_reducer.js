import {
  FETCH_WORKFLOWS,
  FETCH_WORKFLOW,
  SUBMIT_INPUTS
} from '../actions/archimedes_actions';

const archimedesReducer = (archimedes, action) => {
  if (!archimedes) {
    archimedes = {
      workflows: {}
    };
  }
  switch (action.type) {
    case FETCH_WORKFLOWS:
      return {
        workflows: action.workflows ? Object.assign({}, action.workflows) : {}
      };
    case FETCH_WORKFLOW:
      return Object.assign({}, archimedes, {
        workflow: action.workflow
      });
    case SUBMIT_INPUTS:
      // Update the current workflow with status returned from input submission;
      const updatedWorkflow = Object.assign({}, archimedes.workflow);
      Object.keys(action.workflow.steps).forEach((stepName) => {
        updatedWorkflow.steps[stepName] = Object.assign(
          {},
          updatedWorkflow.steps[stepName],
          action.workflow.steps[stepName]
        );
      });
      return Object.assign({}, archimedes, {
        workflow: updatedWorkflow
      });
    default:
      return archimedes;
  }
};

export default archimedesReducer;
