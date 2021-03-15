require('etna-js/spec/setup');

import ReactModal from 'react-modal';
ReactModal.setAppElement('*'); // suppresses modal-related test warnings.

global.ROUTES = {
  workflow: () => `/workflow`,
  fetch_workflows: () => `/api/workflows`,
  submit: (project_name, workflow_name) =>
    `/api/${project_name}/session/${workflow_name}`
};

// Used by plotly.js
global.URL.createObjectURL = jest.fn();
