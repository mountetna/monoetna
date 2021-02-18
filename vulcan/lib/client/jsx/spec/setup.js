require('etna-js/spec/setup');

import ReactModal from 'react-modal';
ReactModal.setAppElement('*'); // suppresses modal-related test warnings.

global.ROUTES = {
  workflow: () => `/${CONFIG.project_name}/workflow`,
  fetch_workflows: () => `/api/${CONFIG.project_name}/workflows`,
  fetch_workflow: (workflow_name, json) =>
    `/api/${CONFIG.project_name}/workflows/${workflow_name}/steps${
      json ? '?format=json' : ''
    }`,
  submit_inputs: (workflow_name) =>
    `/api/${CONFIG.project_name}/workflows/${workflow_name}`
};
