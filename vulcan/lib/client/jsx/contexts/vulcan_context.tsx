import React, {createContext, useRef, useState, useCallback} from 'react';
import VulcanReducer, {defaultVulcanState, VulcanState} from '../reducers/vulcan_reducer';
import {VulcanAction} from '../actions/vulcan_actions';
import {defaultSessionStorageHelpers, useLocalSessionStorage} from './session_storage';
import {defaultApiHelpers, useApi} from './api';
import {defaultInputStateManagement, useInputStateManagement} from './input_state_management';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {defaultConfirmationHelpers, useConfirmation} from './confirmation';
import { useWorkspacesWorkflowLoading } from './workspace_worksflow_loading';

export const defaultContext = {
  state: defaultVulcanState as VulcanState,
  stateRef: {current: defaultVulcanState}, // This would be set with a context dispatch when the provide is actually
                                           // installed.
  dispatch: (a: VulcanAction) => console.warn('action dispatched but not handled', a),
  useActionInvoker: (() => () => null) as typeof useActionInvoker, ...defaultConfirmationHelpers, ...defaultSessionStorageHelpers, ...defaultApiHelpers, ...defaultInputStateManagement,
};

export type VulcanContextData = typeof defaultContext;
export const VulcanContext = createContext(defaultContext);
export type VulcanContext = typeof VulcanContext;
export type ProviderProps = {
  params?: {},
  storage?: typeof localStorage,
  logActions?: boolean,
  wrapper?: [Function, Function],
  children: any
};

// Existing mostly to down type the return result into A from potentially inferring A & B
function withOverrides<A, B extends Partial<A>>(base: A, overrides: B): A {
  const result: A = {...base};
  for (let k in overrides) {
    const v = overrides[k];
    if (v !== undefined) { // @ts-ignore
      result[k] = v;
    }
  }

  return result;
}

export const VulcanProvider = (props: ProviderProps & Partial<VulcanContextData>) => {
  const [start, end] = props.wrapper || [() => null, () => null];
  start();

  try {
    const actionInvokerHelpers = withOverrides({useActionInvoker}, props);
    const invoker = actionInvokerHelpers.useActionInvoker();
    const stateRef = useRef(props.state || defaultContext.state);
    const [state, setState] = useState(stateRef.current);

    const dispatch = useCallback((action: VulcanAction) => {
      if (props.logActions) {
        console.log(action.type, action);
      }

      const result = VulcanReducer(stateRef.current, action);
      stateRef.current = result;
      setState(result);

      return result;
    }, [props.logActions]);

    // Inject data coming from the URL to state.
    if (!!props.params && 'project_name' in props.params && state.projectName!=props.params.project_name) {
      const newResult = {...stateRef.current, projectName: props.params.project_name as string};
      stateRef.current = newResult;
      setState(newResult);
    }
    if (!!props.params && 'workspace_id' in props.params && state.workspaceId!=props.params.workspace_id) {
      const newResult = {...stateRef.current, workspaceId: Number(props.params.workspace_id) as number};
      stateRef.current = newResult;
      setState(newResult);
    }

    const localSessionHelpers = withOverrides(useLocalSessionStorage(state, props), props);
    const apiHelpers = withOverrides(useApi(invoker), props);
    const {showErrors, getWorkflows, getWorkspaces, requestRun, pullRunStatus, getIsRunning} = apiHelpers;
    useWorkspacesWorkflowLoading(state.update_workflows, dispatch, getWorkflows, getWorkspaces, showErrors, state.projectName);
    const confirmationHelpers = withOverrides(useConfirmation(), props);
    const inputHelpers = useInputStateManagement(invoker, dispatch, stateRef);

    return (<VulcanContext.Provider value={{
      state,
      stateRef,
      dispatch,
      ...confirmationHelpers,
      ...actionInvokerHelpers,
      ...localSessionHelpers,
      ...apiHelpers,
      ...inputHelpers,
    }}>
      {props.children}
    </VulcanContext.Provider>);
  } finally {
    end();
  }
};
