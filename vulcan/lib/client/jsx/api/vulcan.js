import * as _ from 'lodash';

import {
  checkStatus,
  handleFetchSuccess,
  handleFetchError,
  headers,
  isJSON
} from 'etna-js/utils/fetch';
import {shouldDownloadStepData} from '../utils/workflow';

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

const rawVulcanGet = (endpoint) => {
  return fetch(endpoint, {
    method: 'GET',
    credentials: 'include',
    headers: headers('json')
  });
};

const vulcanGet = (endpoint) => {
  return rawVulcanGet(endpoint).then(checkStatus);
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

export const submit = (context) => {
  let {
    workflow,
    session,
    status,
    pathIndex,
    setSession,
    setStatus,
    setData
  } = context;

  let dataUrls = [];

  // Use a deep clone because we need to check if the download URL
  //   changed, later.
  let oldStatus = _.cloneDeep(status);

  return vulcanPost(vulcanPath(ROUTES.submit(workflow.name)), {
    inputs: session.inputs,
    key: session.key
  }).then(handleFetchSuccess)
    .then((response) => {
      setSession(response.session);
      setStatus(response.status);

      // Calculate this locally because useContext
      //   updates async?
      // Right -- we should move away from reducers and immutable state, and instead just have service objects
      // that are mutable  being shared via the context, that way once we have a reference to that mutable
      // service object, state changes can be shared via that shared reference, rather than being forced to await the
      // new context value being propagated up.
      let updatedStatus = response.status.map((newPath, newPathIndex) => {
        return newPath.map((newStep, newStepIndex) => {
          let oldStep = {};
          if (
            oldStatus &&
            oldStatus[newPathIndex] &&
            oldStatus[newPathIndex][newStepIndex]
          )
            oldStep = oldStatus[newPathIndex][newStepIndex];

          return {
            ...oldStep,
            ...newStep
          };
        });
      });

      let dataRequests = [];
      updatedStatus[pathIndex].forEach((step, stepIndex) => {
        if (step.downloads) {
          Object.keys(step.downloads).forEach((download) => {
            if (
              downloadUrlUpdated(
                oldStatus[pathIndex][stepIndex],
                step,
                download
              ) &&
              shouldDownloadStepData({workflow, pathIndex, stepIndex})
            ) {
              let dataUrl = step.downloads[download];
              dataUrls.push(dataUrl);
              dataRequests.push(
                // Don't want checkStatus here so that we
                //  can call Promise.all(), just the raw fetch.
                rawVulcanGet(dataUrl)
              );
            }
          });
        }
      });
      return Promise.all(dataRequests);
    })
    .then((downloads) => {
      let extractData = [];

      downloads.forEach((response) => {
        let data = isJSON(response) ? response.json() : response.text();

        extractData.push(data);
      });

      return Promise.all(extractData);
    })
    .then((extractions) => {
      extractions.forEach((datum, index) => {
        // We don't currently have accurate Content-Type headers
        //   from the response, so we'll try here if it's a
        //   JSON string or not.
        try {
          datum = JSON.parse(datum);
        } catch (e) {}

        setData(dataUrls[index], datum);
      });
      return Promise.resolve();
    })
    .catch(handleFetchError);
};

// export for testing
export const downloadUrlUpdated = (oldStep, newStep, downloadKey) => {
  return (
    !(newStep.data && newStep.data[downloadKey]) ||
    newStep.downloads[downloadKey] !== oldStep.downloads[downloadKey]
  );
};

export const getData = (url) => {
  return vulcanGet(url).then(handleFetchSuccess).catch(handleFetchError);
};
