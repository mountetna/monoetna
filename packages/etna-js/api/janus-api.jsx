export const getProjects = () => {
  let projects = require('./projects.json');

  return new Promise((resolve, reject) => {
    resolve(projects);
  });
};
