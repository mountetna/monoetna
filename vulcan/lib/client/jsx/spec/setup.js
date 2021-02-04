require('etna-js/spec/setup');

import ReactModal from 'react-modal';
ReactModal.setAppElement('*'); // suppresses modal-related test warnings.

global.ROUTES = {
  workflow: () => `/${CONFIG.project_name}/workflow`,
  fetch_workflows: () => '/api/workflows',
  fetch_workflow: (workflow_name) => `/api/workflows/${workflow_name}/steps`,
  fetch_pools: (workflow_name) => `/api/workflows/${workflow_name}/pools`,
  submit_inputs: (workflow_name) => `/api/workflows/${workflow_name}`
};
