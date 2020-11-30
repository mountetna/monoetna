import {FETCH_PROJECTS} from '../actions/janus-actions';

const janusReducer = (janus, action) => {
  if (!janus) {
    janus = {
      projects: []
    };
  }
  switch (action.type) {
    case FETCH_PROJECTS:
      return {
        projects: [...action.projects]
      };
    default:
      return janus;
  }
};

export default janusReducer;
