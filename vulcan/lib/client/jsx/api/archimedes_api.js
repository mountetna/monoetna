import {checkStatus, headers} from 'etna-js/utils/fetch';

import CwlSerializer from '../serializers/cwl';

const archimedesPath = (endpoint) => `${CONFIG.archimedes_host}/${endpoint}`;

const archimedesPost = (endpoint, exchange, params) => {
  return exchange
    .fetch(archimedesPath(endpoint), {
      method: 'POST',
      credentials: 'include',
      headers: headers('json'),
      body: JSON.stringify({
        ...params,
        project_name: CONFIG.project_name
      })
    })
    .then(checkStatus);
};

const archimedesGet = (endpoint, exchange, params) => {
  return exchange
    .fetch(archimedesPath(endpoint), {
      method: 'GET',
      credentials: 'include',
      headers: headers('json')
    })
    .then(checkStatus);
};

export const getWorkflows = (exchange) => {
  // TODO: update this per real Archimedes endpoint
  return archimedesGet('api/workflows', exchange);
};

export const getWorkflow = (workflow_name, exchange) => {
  // TODO: update this per real Archimedes endpoint
  return archimedesGet(`api/workflows/${workflow_name}`, exchange).then(
    (response) => {
      return new Promise((resolve, reject) => {
        resolve(new CwlSerializer(response));
      });
    }
  );
};

export const submitInputs = (workflow_name, step, inputs, exchange) => {
  // TODO: remove the "step" stub for real Archimedes
  return archimedesPost(
    `api/workflows/${workflow_name}/${step}`,
    exchange,
    inputs
  );
};

export const getData = (workflow_name, data, exchange) => {
  // TODO: update this per real Archimedes endpoint
  return archimedesGet(`api/workflows/${workflow_name}/${data}`, exchange);
};
