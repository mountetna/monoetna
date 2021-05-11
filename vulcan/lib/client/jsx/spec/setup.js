require('etna-js/spec/setup');

import ReactModal from 'react-modal';
ReactModal.setAppElement('*'); // suppresses modal-related test warnings.

import 'jest-canvas-mock';

global.ROUTES = {
  workflow: () => `/workflow`,
  workflow_vignette: (workflow_name) => `/workflow/${workflow_name}/vignette`,
  fetch_workflows: () => `/api/workflows`,
  submit: (project_name, workflow_name) =>
    `/api/${project_name}/session/${workflow_name}`,
  status: (project_name, workflow_name) =>
    `/api/${project_name}/session/${workflow_name}/status`,
};

// Used by plotly.js
global.URL.createObjectURL = jest.fn();
