import * as _ from 'lodash';

import {
  checkStatus,
  handleFetchSuccess,
  handleFetchError,
  headers,
  isJSON
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
  })
    .then(handleFetchSuccess)
    .then((response) => {
      setSession(response.session);
      setStatus(response.status);
      // Fetch data and update Context

      // Calculate this locally because useContext
      //   updates async?
      let updatedStatus = [...oldStatus].map((oldPath, oldPathIndex) => {
        return oldPath.map((oldStep, oldStepIndex) => {
          return {...oldStep, ...response.status[oldPathIndex][oldStepIndex]};
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
              )
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
