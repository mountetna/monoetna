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
      return archimedes;
    default:
      return archimedes;
  }
};

export default archimedesReducer;
