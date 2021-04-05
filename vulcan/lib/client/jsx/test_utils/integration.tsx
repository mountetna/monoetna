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
import {patchInputs, setDownloadedData, setStatus, VulcanAction} from "../actions/vulcan";
import {splitSource, statusOfStep, stepOfSource} from "../selectors/workflow_selectors";
import {defaultStepStatus, StepStatus} from "../api_types";
import {createStatusFixture, createStepStatusFixture} from "./fixtures";

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

  return {node, updateMatching, contextData, reduxState, replaceOverrides, dispatch, setData};

  function setData(sourceName: string, value: any) {
    const workflow = contextData.state.workflow;
    if (!workflow) throw new Error('Workflow must be set!');

    let url: string | undefined;

    const [stepName, outputName] = splitSource(sourceName);
    if (!stepName) throw new Error('Cannot setData for primary inputs');

    let status = statusOfStep(stepName, contextData.state.status) || defaultStepStatus;
    const {downloads} = status;
    if (downloads) {
      url = downloads[splitSource(sourceName)[1]]
    }
    if (!url) url = "https://" + sourceName;


    contextData.dispatch(setStatus(createStatusFixture(workflow, createStepStatusFixture({
      ...status,
      name: stepName,
      status: 'complete',
      downloads: {...status.downloads, [outputName]: url}
    }))));

    contextData.dispatch(setDownloadedData(url, value));
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
    await act(async function () {
      contextData.dispatch(action);
    })
  }
}