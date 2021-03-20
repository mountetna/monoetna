import React, {useReducer, createContext} from 'react';
import {
  SET_DATA,
  SET_WORKFLOW,
  SET_WORKFLOWS,
  SET_STATUS,
  SET_STEP,
  SET_PATH,
  SET_SESSION,
  SET_INPUTS,
  SET_CALCULATING
} from '../actions/vulcan';
import VulcanReducer from '../reducers/vulcan';
import {allInputsDefined, localStorageKey} from '../utils/workflow';

export const VulcanContext = createContext();

export const VulcanProvider = (props) => {
  const initialState = {
    workflows: {},
    workflow: {},
    pathIndex: 0,
    stepIndex: null,
    session: {inputs: {}},
    status: [[]],
    calculating: false,
    ...props.state // useful for testing
  };

  const [state, dispatch] = useReducer(VulcanReducer, initialState);

  const setWorkflows = (workflows) => {
    dispatch({type: SET_WORKFLOWS, workflows});
  };

  const setWorkflow = (workflow) => {
    dispatch({type: SET_WORKFLOW, workflow});
  };

  const setStatus = (status) => {
    dispatch({type: SET_STATUS, status});
  };

  const setData = (url, data) => {
    dispatch({type: SET_DATA, url, data});
  };

  const setPathIndex = (pathIndex) => {
    dispatch({type: SET_PATH, pathIndex});
  };

  const setStepIndex = (stepIndex) => {
    dispatch({type: SET_STEP, stepIndex});
  };

  const setSession = (session) => {
    if (state.workflow && state.workflow.name)
      localStorage.setItem(
        localStorageKey(state.workflow),
        JSON.stringify(session)
      );
    dispatch({type: SET_SESSION, session});
  };

  const setInputs = (inputs) => {
    dispatch({type: SET_INPUTS, inputs});
  };

  const setCalculating = (calculating) => {
    dispatch({type: SET_CALCULATING, calculating});
  };

  const getLocalSession = (workflow) => {
    let storedSession = localStorage.getItem(localStorageKey(workflow));
    if (!storedSession) return Promise.resolve(null);

    // We now need to check if the input names have changed.
    // If all the workflow's primary inputs are NOT present
    //   in the stored session, we'll return `null` and get
    //   a new session.
    storedSession = JSON.parse(storedSession);
    if (!allInputsDefined(workflow, storedSession.inputs)) {
      localStorage.removeItem(localStorageKey(workflow));
      return Promise.resolve(null);
    }

    return Promise.resolve(storedSession);
  };

  return (
    <VulcanContext.Provider
      value={{
        workflows: state.workflows,
        workflow: state.workflow,
        pathIndex: state.pathIndex,
        stepIndex: state.stepIndex,
        session: state.session,
        status: state.status,
        calculating: state.calculating,
        setWorkflows,
        setWorkflow,
        setData,
        setStatus,
        setPathIndex,
        setStepIndex,
        setSession,
        setInputs,
        setCalculating,
        getLocalSession
      }}
    >
      {props.children}
    </VulcanContext.Provider>
  );
};
