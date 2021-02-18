import {
  checkStatus,
  handleFetchSuccess,
  handleFetchError,
  headers
} from 'etna-js/utils/fetch';

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
  return vulcanGet(ROUTES.fetch_workflows())
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

export const submit = (workflow_name, inputs) => {
  return vulcanPost(ROUTES.submit_inputs(workflow_name), inputs)
    .then(handleFetchSuccess)
    .catch(handleFetchError);
};

export const getData = (url) => {
  return vulcanGet(url).then(handleFetchSuccess).catch(handleFetchError);
};
