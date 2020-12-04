import {getProjects} from '../api/janus-api';

export const FETCH_PROJECTS = 'FETCH_PROJECTS';

export const fetchProjectsAction = () => {
  return (dispatch) => {
    getProjects().then((projects) => {
      dispatch({
        type: FETCH_PROJECTS,
        projects: projects
      });
    });
  };
};
