import React, {createContext, useRef, useState, useCallback} from 'react';
import VulcanReducer, {defaultVulcanState, VulcanState} from '../reducers/vulcan_reducer';
import {VulcanAction} from '../actions/vulcan_actions';
import {defaultSessionStorageHelpers, useLocalSessionStorage} from './session_storage';
import {getWorkflows, getWorkspaces} from './api';
import {defaultInputStateManagement, useInputStateManagement} from './input_state_management';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {defaultConfirmationHelpers, useConfirmation} from './confirmation';
import { useWorkspacesWorkflowLoading } from './workspace_worksflow_loading';

export const defaultContext = {
  state: defaultVulcanState as VulcanState,
};

export type VulcanContextData = typeof defaultContext;
export const VulcanContext = createContext(defaultContext);
export type VulcanContext = typeof VulcanContext;
export type ProviderProps = {
  params?: {},
  children: any
};

export const VulcanProvider = (props: ProviderProps & Partial<VulcanContextData>) => {
  const [ workflows, setWorkflows ] = useState([]);
  const [ workspaces, setWorkspaces ] = useState([]);

  const projectName = props.params.project_name;

  React.useEffect(
    () => {
      console.log("Getting "+projectName);
      getWorkspaces(projectName).then(({workspaces}) => setWorkspaces(workspaces));
      getWorkflows(projectName).then(({workflows}) => setWorkflows(workflows));
    }, [projectName]
  );

  return (<VulcanContext.Provider value={{
    state: {
      projectName,
      workflows,
      workspaces
    }
  }}>
    {props.children}
  </VulcanContext.Provider>);
};
