import {json_get} from '../utils/fetch';

export const getProjects = () => {
  return json_get(`${CONFIG.janus_host}/api/user/projects`).then((projects) => {
    return new Promise((resolve, reject) => {
      resolve(projects.projects);
    });
  });
};
