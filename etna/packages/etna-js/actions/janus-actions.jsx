import {getProjects} from '../api/janus-api';
import {handleFetchError} from '../utils/fetch';

export const FETCH_PROJECTS = 'FETCH_PROJECTS';

export const fetchProjectsAction = () => {
  return (dispatch) => {
    getProjects().then((projects) => {
      dispatch({
        type: FETCH_PROJECTS,
        projects: projects
      });
    }).catch(handleFetchError);
  };
};
