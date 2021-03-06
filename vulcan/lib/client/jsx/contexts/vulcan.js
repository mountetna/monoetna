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
    dispatch({type: SET_SESSION, session});
  };

  const setInputs = (inputs) => {
    dispatch({type: SET_INPUTS, inputs});
  };

  const setCalculating = (calculating) => {
    dispatch({type: SET_CALCULATING, calculating});
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
        setCalculating
      }}
    >
      {props.children}
    </VulcanContext.Provider>
  );
};
