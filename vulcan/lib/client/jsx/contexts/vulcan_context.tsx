import React, {useReducer, createContext, useRef} from 'react';
import VulcanReducer, {defaultVulcanState, VulcanState} from '../reducers/vulcan_reducer';
import {VulcanAction} from "../actions/vulcan_actions";
import {defaultSessionStorageHelpers, useLocalSessionStorage} from "./session_storage";
import {defaultApiHelpers, useApi} from "./api";
import {useWorkflowsLoading} from "./workflows_loading";
import {useDataBuffering} from "./data_buffering";
import {defaultSessionSyncHelpers, useSessionSync} from "./session_sync";
import {defaultInputStateManagement, useInputStateManagement} from "./input_state_management";
import {useActionInvoker} from "etna-js/hooks/useActionInvoker";
import {defaultConfirmationHelpers, useConfirmation} from "./confirmation";

export const defaultContext = {
  state: defaultVulcanState as VulcanState,
  stateRef: {current: defaultVulcanState}, // This would be set with a context dispatch when the provide is actually
                                           // installed.
  dispatch: (a: VulcanAction) => console.warn('action dispatched but not handled', a),
  useActionInvoker: (() => () => null) as typeof useActionInvoker, ...defaultConfirmationHelpers, ...defaultSessionStorageHelpers, ...defaultApiHelpers, ...defaultSessionSyncHelpers, ...defaultInputStateManagement,
}

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
    const [state, dispatch] = useReducer(function (state: VulcanState, action: VulcanAction) {
      if (props.logActions) {
        console.log(action.type, action);
      }

      const result = VulcanReducer(state, action);

      stateRef.current = result;
      return result;
    }, stateRef.current);
    const localSessionHelpers = withOverrides(useLocalSessionStorage(state, props), props);
    const apiHelpers = withOverrides(useApi(invoker), props);
    const {showErrors, getData, getWorkflows, pollStatus, postInputs} = apiHelpers;
    const sessionSyncHelpers = withOverrides(
      useSessionSync(stateRef, showErrors, pollStatus, postInputs, dispatch),
      props
    );
    useDataBuffering(state, dispatch, showErrors, getData);
    useWorkflowsLoading(JSON.stringify(props.params), dispatch, getWorkflows, showErrors);
    const confirmationHelpers = withOverrides(useConfirmation(), props);
    const inputHelpers = useInputStateManagement(invoker, dispatch, sessionSyncHelpers.requestPoll, stateRef);

    return (<VulcanContext.Provider value={{
      state,
      stateRef,
      dispatch, ...confirmationHelpers, ...actionInvokerHelpers, ...localSessionHelpers, ...apiHelpers, ...sessionSyncHelpers, ...inputHelpers,
    }}>
      {props.children}
    </VulcanContext.Provider>);
  } finally {
    end();
  }
};
