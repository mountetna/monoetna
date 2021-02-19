import {
  checkStatus,
  handleFetchSuccess,
  handleFetchError,
  headers
} from 'etna-js/utils/fetch';

const vulcanPath = (endpoint) => `${CONFIG.vulcan_host}${endpoint}`;

const vulcanPost = (endpoint, params) => {
  return fetch(endpoint, {
    method: 'POST',
    credentials: 'include',
    headers: headers('json'),
    body: JSON.stringify({
      ...params
    })
  }).then(checkStatus);
};

const vulcanGet = (endpoint) => {
  return fetch(endpoint, {
    method: 'GET',
    credentials: 'include',
    headers: headers('json')
  }).then(checkStatus);
};

export const getWorkflows = () => {
  return vulcanGet(vulcanPath(ROUTES.fetch_workflows()))
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};

// export const getWorkflow = (workflow_name, json = false) => {
//   // NOTE: response might be a YAML doc, for editing purposes.
//   // TODO: update this per real Vulcan endpoint
//   return vulcanGet(ROUTES.fetch_workflow(workflow_name, json))
//     .then(handleFetchSuccess)
//     .catch(handleFetchError);
// };

export const submit = (workflow_name, inputs, key) => {
  return vulcanPost(vulcanPath(ROUTES.submit(workflow_name)), {inputs, key})
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};

export const getSession = (workflow_name) => {
  // A "blank" POST to submit generates and returns the
  //   session, which we'll need for subsequent input
  //   submit actions.
  return vulcanPost(vulcanPath(ROUTES.submit(workflow_name)))
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};

export const getData = (url) => {
  return vulcanGet(url).then(handleFetchSuccess).catch(handleFetchError);
};
