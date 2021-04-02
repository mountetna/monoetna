import {ReactElement, useContext} from "react";
import {act, create} from "react-test-renderer";
import {
  defaultContext,
  ProviderProps,
  VulcanContext,
  VulcanContextData,
  VulcanProvider
} from "../contexts/vulcan_context";
import {useReduxState} from "etna-js/hooks/useReduxState";
import {Provider} from "react-redux";
import {VulcanStore} from "../vulcan_store";
import * as React from "react";
import {Store} from "redux";
import {setInputs, VulcanAction} from "../actions/vulcan";

function injectContextAgent(Element: () => ReactElement | null, overrides: Partial<ProviderProps & VulcanContextData>) {
  return <VulcanProvider {...overrides}>
    <Element/>
  </VulcanProvider>
}

export function integrateElement(Element: () => ReactElement | null, {
  store = VulcanStore(),
  wrapper = injectContextAgent,
  providerOverrides = {},
}: { store?: Store, wrapper?: typeof injectContextAgent, providerOverrides?: Partial<ProviderProps & VulcanContextData> } = {}) {
  let contextData: VulcanContextData = defaultContext;
  let reduxState: any = {};
  let waiters: Function[] = [];

  function updateMatching(pred: () => boolean): Promise<void> {
    return new Promise((resolve, reject) => {
      let timedOut = false;

      const t = setTimeout(function () {
        timedOut = true;
        reject('timeout');
      }, 3000);

      waiters.push(function () {
        if (!pred()) return false;
        if (timedOut) return true;

        clearTimeout(t);
        resolve();
        return true;
      });
    })
  }

  const node = create(
      <Provider store={store}>
        {wrapper(TestComponent, providerOverrides)}
      </Provider>
  )

  return {node, updateMatching, contextData, reduxState, setInput, replaceOverrides, dispatch}

  // Utility functions
  function setInput(inputName: string, value: any) {
    contextData.dispatch(setInputs({...contextData.state.inputs, [inputName]: value}));
  }

  function replaceOverrides(overrides: Partial<ProviderProps & VulcanContextData>) {
    node.update(<Provider store={store}>
      {wrapper(TestComponent, overrides)}
    </Provider>)
  }

  function TestComponent() {
    // Keep the context data inside the shared reference returned back to the called.
    // Because context data's shape never changes, it is safe to simply call assign without
    // clearing it first.
    Object.assign(contextData, useContext(VulcanContext));

    // Ensure no polution of the redux state -- the shape of which can change between calls.
    Object.keys(reduxState).forEach(k => delete reduxState[k]);
    Object.assign(reduxState, useReduxState());

    // Check each pending waiter with the update received on the component, and if it has completed, remove it.
    // We go in reverse order to prevent the issues involved with deleting during iteration.
    for (let i = waiters.length - 1; i >= 0; --i) {
      if (waiters[i]()) waiters.splice(i, 1);
    }

    return <Element/>;
  }

  async function dispatch(action: VulcanAction) {
    await act(async function() {
      contextData.dispatch(action);
    })
  }
}