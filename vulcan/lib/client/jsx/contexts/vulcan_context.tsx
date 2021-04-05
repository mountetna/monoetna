import React, {useReducer, createContext, useRef, useEffect} from 'react';
import VulcanReducer, {defaultVulcanState, VulcanState} from '../reducers/vulcan_reducer';
import {VulcanAction} from "../actions/vulcan";
import {defaultSessionStorageHelpers, useLocalSessionStorage} from "./session_storage";
import {defaultApiHelpers, useApi} from "./api";
import {useWorkflowsLoading} from "./workflows_loading";
import {defaultDataBufferingHelpers, useDataBuffering} from "./data_buffering";
import {defaultSessionSyncHelpers, useSessionSync} from "./session_sync";
import {useInputStateManagement} from "./input_state_management";
import {useActionInvoker} from "etna-js/hooks/useActionInvoker";

export const defaultContext = {
  state: defaultVulcanState as VulcanState,
  stateRef: {current: defaultVulcanState},
  // This would be set with a context dispatch when the provide is actually installed.
  dispatch: (a: VulcanAction) => console.warn('action dispatched but not handled', a),
  useActionInvoker: (() => () => null) as typeof useActionInvoker,
  ...defaultSessionStorageHelpers,
  ...defaultApiHelpers,
  ...defaultSessionSyncHelpers,
  ...defaultDataBufferingHelpers,
}

export type VulcanContextData = typeof defaultContext;
export const VulcanContext = createContext(defaultContext);
export type VulcanContext = typeof VulcanContext;
export type ProviderProps = { params?: {}, storage?: typeof localStorage, children: any };

// Existing mostly to down type the return result into A from potentially inferring A & B
function withOverrides<A, B extends Partial<A>, K extends keyof A>(base: A, overrides: B): A {
  return {...base, ...overrides};
}

export const VulcanProvider = (props: ProviderProps & Partial<VulcanContextData>) => {
  const actionInvokerHelpers = withOverrides({useActionInvoker}, props);
  const stateRef = useRef(props.state || defaultContext.state);
  const [state, dispatch] = useReducer(function (state: VulcanState, action: VulcanAction) {
    const result = VulcanReducer(state, action);
    stateRef.current = result;
    return result;
  }, stateRef.current);
  const localSessionHelpers = withOverrides(useLocalSessionStorage(state, props), props);
  const apiHelpers = withOverrides(useApi(actionInvokerHelpers.useActionInvoker()), props);
  const {scheduleWork, getData, getWorkflows, pollStatus, postInputs} = apiHelpers;
  const sessionSyncHelpers = withOverrides(useSessionSync(stateRef, scheduleWork, pollStatus, postInputs, dispatch), props);
  const dataBufferingHelpers = useDataBuffering(state, dispatch, scheduleWork, getData);
  useWorkflowsLoading(JSON.stringify(props.params), dispatch, getWorkflows, scheduleWork);
  useInputStateManagement(state, dispatch, sessionSyncHelpers.statusIsFresh);


  return (
      <VulcanContext.Provider value={{
        state,
        stateRef,
        dispatch,
        ...actionInvokerHelpers,
        ...localSessionHelpers,
        ...apiHelpers,
        ...sessionSyncHelpers,
        ...dataBufferingHelpers,
      }}>
        {props.children}
      </VulcanContext.Provider>
  );
};
