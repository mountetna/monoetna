import {Exchange} from 'etna-js/actions/exchange_actions';
import {showMessages} from 'etna-js/actions/message_actions';
import {handleFetchError, checkFetchSuccess} from 'etna-js/utils/fetch';
import {
  getData,
  getWorkflow,
  getWorkflows,
  submitInputs
} from '../api/archimedes_api';
import CwlSerializer from '../serializers/cwl';

export const FETCH_WORKFLOWS = 'FETCH_WORKFLOWS';
export const FETCH_WORKFLOW = 'FETCH_WORKFLOW';
export const SUBMIT_INPUTS = 'SUBMIT_INPUTS';
export const FETCH_DATA = 'FETCH_DATA';

const fetchWorkflowsAction = (workflows) => {
  return {
    type: FETCH_WORKFLOWS,
    workflows
  };
};

const fetchWorkflowAction = (workflow) => {
  return {
    type: FETCH_WORKFLOW,
    workflow
  };
};

const submitInputAction = (status) => {
  return {
    type: SUBMIT_INPUTS,
    status
  };
};

const fetchDataAction = (url) => {
  return {
    type: FETCH_DATA,
    url
  };
};

export const fetchWorkflows = (dispatch) => {
  const exchange = new Exchange(dispatch, 'archimedes-fetch-workflows');
  return getWorkflows(exchange)
    .then((response) => {
      return checkFetchSuccess(dispatch, response);
    })
    .then((response) => {
      dispatch(fetchWorkflowsAction(response.workflows));
    })
    .catch((e) => {
      handleFetchError(dispatch, e);
    });
};

export const fetchWorkflow = (workflow_name) => {
  return (dispatch) => {
    const exchange = new Exchange(dispatch, 'archimedes-fetch-workflow');
    return getWorkflow(workflow_name, exchange)
      .then((response) => {
        return checkFetchSuccess(dispatch, response);
      })
      .then((response) => {
        // response should be a JSON doc here by now
        dispatch(fetchWorkflowAction(response));
      })
      .catch((e) => {
        handleFetchError(dispatch, e);
      });
  };
};

export const submitWorkflowInputs = (inputs) => {
  return (dispatch) => {
    const exchange = new Exchange(dispatch, 'archimedes-submit-inputs');
    return submitInputs(workflow_name, inputs, exchange)
      .then((response) => {
        return checkFetchSuccess(dispatch, response);
      })
      .then((response) => {
        // response should be the current status, given the set of inputs?
        dispatch(submitInputAction(response));
      })
      .catch((e) => {
        handleFetchError(dispatch, e);
      });
  };
};
