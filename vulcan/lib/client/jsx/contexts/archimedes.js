import React, {useReducer, createContext} from 'react';
import {
  SET_DATA,
  SET_WORKFLOW,
  SET_WORKFLOWS,
  SET_STATUS
} from '../actions/archimedes';
import ArchimedesReducer from '../reducers/archimedes';

export const ArchimedesContext = createContext();

export const ArchimedesProvider = (props) => {
  const initialState = {
    workflows: {},
    workflow: {}
  };

  const [state, dispatch] = useReducer(ArchimedesReducer, initialState);

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

  return (
    <ArchimedesContext.Provider
      value={{
        workflows: state.workflows,
        workflow: state.workflow,
        setWorkflows,
        setWorkflow,
        setData,
        setStatus
      }}
    >
      {props.children}
    </ArchimedesContext.Provider>
  );
};
