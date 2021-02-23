require('etna-js/spec/setup');

import ReactModal from 'react-modal';
ReactModal.setAppElement('*'); // suppresses modal-related test warnings.

global.ROUTES = {
  workflow: () => `/${CONFIG.project_name}/workflow`,
  fetch_workflows: () => `/api/${CONFIG.project_name}/workflows`,
  submit: (workflow_name) =>
    `/api/${CONFIG.project_name}/session/${workflow_name}`
};
