import 'regenerator-runtime/runtime';
// suppresses a specific React warning (comment raf out to bring it back)
import nock from 'nock';
import {stringifyRequest} from 'nock/lib/common';
const raf = (global.requestAnimationFrame = (cb) => setTimeout(cb, 0));

jest.mock('webworkify-webpack');

const Enzyme = require('enzyme');
const EnzymeAdapter = require('enzyme-adapter-react-16');
// Setup enzyme's react adapter
Enzyme.configure({adapter: new EnzymeAdapter()});

global.fetch = require('isomorphic-fetch');

nock.disableNetConnect();
nock.emitter.on('no match', (...args) => {
  // nock has inconsistent parameters for this event based on how the match failed, which results in us having complex
  // unpacking logic.
  let [_, options, bodyString] = args;
  if (options == null) {
    ({options} = args[0]);
    if (options == null) {
      options = args[0];
    }
  }

  const reqString = stringifyRequest(options, bodyString);
  console.error('got request with no match', reqString);
});

// nodejs equivalent
global.FormData = URLSearchParams;

global.CONFIG = {
  baseURL: 'http://localhost'
};

// For plots/actions/__tests__/manifest.test.js
global.Routes = {
  fetch_manifests_path: (projectName) =>
    `http://localhost/${projectName}/manifests`,
  destroy_manifest_path: (projectName, manifestId) =>
    `http://localhost/${projectName}/manifests/destroy/${manifestId}`,
  create_manifest_path: (projectName) =>
    `http://localhost/${projectName}/manifests/create`,
  update_manifest_path: (projectName, manifestId) =>
    `http://localhost/${projectName}/manifests/update/${manifestId}`,
  plots_fetch_path: (project_name) => `http://localhost/${project_name}/plots`,
  view_path: (project_name, model_name) =>
    `http://localhost/${project_name}/view/${model_name}`,
  browse_model_path: (project_name, model_name, record_name) =>
    `http://localhost/${project_name}/browse/${model_name}/${record_name}`
};
