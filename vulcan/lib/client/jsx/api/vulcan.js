import {
  checkStatus,
  handleFetchSuccess,
  handleFetchError,
  headers
} from 'etna-js/utils/fetch';

const vulcanPath = (endpoint) => `${CONFIG.vulcan_host}${endpoint}`;

const vulcanPost = (endpoint, params) => {
  return fetch(vulcanPath(endpoint), {
    method: 'POST',
    credentials: 'include',
    headers: headers('json'),
    body: JSON.stringify({
      ...params
    })
  }).then(checkStatus);
};

const vulcanGet = (endpoint) => {
  return fetch(vulcanPath(endpoint), {
    method: 'GET',
    credentials: 'include',
    headers: headers('json')
  }).then(checkStatus);
};

export const getWorkflows = () => {
  // TODO: update this per real Vulcan endpoint
  return vulcanGet(ROUTES.fetch_workflows())
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};

export const getWorkflow = (workflow_name, json = false) => {
  // NOTE: response might be a YAML doc, for editing purposes.
  // TODO: update this per real Vulcan endpoint
  return vulcanGet(ROUTES.fetch_workflow(workflow_name, json))
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};

export const submitInputs = (workflow_name, inputs) => {
  // TODO: remove the "step" stub for real Vulcan
  return vulcanPost(ROUTES.submit_inputs(workflow_name), inputs)
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};

export const getData = (url) => {
  // TODO: update this per real Vulcan endpoint
  return vulcanGet(url).then(handleFetchSuccess).catch(handleFetchError);
};
