import {
  checkStatus,
  handleFetchSuccess,
  handleFetchError,
  headers
} from 'etna-js/utils/fetch';

import CwlSerializer from '../serializers/cwl';

const vulcanPath = (endpoint) => `${CONFIG.vulcan_host}${endpoint}`;

const vulcanPost = (endpoint, params) => {
  return fetch(vulcanPath(endpoint), {
    method: 'POST',
    credentials: 'include',
    headers: headers('json'),
    body: JSON.stringify({
      ...params,
      project_name: CONFIG.project_name
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

export const getWorkflow = (workflow_name) => {
  // TODO: update this per real Vulcan endpoint
  return vulcanGet(ROUTES.fetch_workflow(workflow_name))
    .then(handleFetchSuccess)
    .then((response) => {
      // response should be a YAML doc
      return Promise.resolve(new CwlSerializer(response).json);
    })
    .catch(handleFetchError);
};

export const submitInputs = (workflow_name, inputs) => {
  // TODO: remove the "step" stub for real Vulcan
  // We'll have to convert inputs into a Hash for the backend,
  //   most likely.
  const inputsHash = inputs.reduce((result, stepInput) => {
    result[stepInput.name] = stepInput;
    return result;
  }, {});

  // NOTE: the returned data from the server may have to be
  //   massaged by the reducer / consumer to fit back into the Array of steps.
  return vulcanPost(ROUTES.submit_inputs(workflow_name), inputsHash)
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};

export const getData = (url) => {
  // TODO: update this per real Vulcan endpoint
  return vulcanGet(url);
};
