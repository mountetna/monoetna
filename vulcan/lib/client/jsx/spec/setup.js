require('etna-js/spec/setup');

import ReactModal from 'react-modal';
ReactModal.setAppElement('*'); // suppresses modal-related test warnings.

global.ROUTES = {
  workflow: (project_name) => `/${project_name}/workflow`,
  fetch_steps: (project_name, workflow_name) =>
    `/api/workflows/${workflow_name}/steps`,
  fetch_pools: (project_name, workflow_name) =>
    `/api/workflows/${workflow_name}/pools`,
  submit_steps: (project_name, workflow_name, status) =>
    `/api/workflows/${workflow_name}/${status}`
};
