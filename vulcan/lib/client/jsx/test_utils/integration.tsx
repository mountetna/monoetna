import {ReactElement, useContext} from "react";
import {act, create} from "react-test-renderer";
import {
  defaultContext, ProviderProps, VulcanContext, VulcanContextData, VulcanProvider
} from "../contexts/vulcan_context";
import {Provider} from "react-redux";
import {VulcanStore} from "../vulcan_store";
import * as React from "react";
import {createFakeStorage} from "./mocks";
import {Cancellable} from "etna-js/utils/cancellable";

export function integrateElement(
  element: ReactElement | null = null,
  providerOverrides: Partial<ProviderProps & VulcanContextData> = {},
) {
  // Default provide mocks should proceed the spread below.
  providerOverrides = {
    getWorkflows: defaultContext.getWorkflows,
    getData: defaultContext.getData,
    pollStatus: defaultContext.pollStatus,
    postInputs: defaultContext.postInputs,
    storage: createFakeStorage(),
    ...providerOverrides,
  }

  const store = VulcanStore();

  let preHook: Function | undefined = undefined;
  let lastPreHook: Function | undefined = preHook;
  let prehookResult = [undefined] as [any];
  let renderIdx = 0;
  let PreHookContainer = function PreHookContainer() {
    if(preHook) prehookResult[0] = preHook();
    return null;
  }

  const node = create(rerender())

  return {node, runHook, update};

  // Regenerates the react element tree to create or update.
  // Manages & injects the prehook as well as the inner element.
  function rerender() {
    if (preHook !== lastPreHook) {
      renderIdx++;
      lastPreHook = preHook;
      prehookResult[0] = undefined;
      PreHookContainer = function PreHookContainer() {
        if(preHook) prehookResult[0] = preHook();
        return null;
      }
    }

    return <Provider store={store}>
      <VulcanProvider {...providerOverrides}>
        <PreHookContainer key={renderIdx + "-prehook"}/>
        {element}
      </VulcanProvider>
    </Provider>;
  }

  function* runHook<T>(hook: () => T): Generator<unknown, T, T> {
    yield function() {
      preHook = hook;
      node.update(rerender());
    }

    return prehookResult[0] as T;
  }

  function* getContext(): Generator<unknown, VulcanContextData, VulcanContextData> {
    return yield* runHook(function() {
      return useContext(VulcanContext);
    })
  }

  function* update<T = unknown>(
    {updateOverrides, replaceElement}: {
      updateOverrides?: (curOverrides: Partial<ProviderProps & VulcanContextData>) =>
        Partial<ProviderProps & VulcanContextData>,
      replaceElement?: ReactElement | null,
    } = {}
  ) {
    if (updateOverrides) {
      providerOverrides = updateOverrides(providerOverrides);
    }

    if (replaceElement !== undefined) {
      element = replaceElement;
    }

    return yield* getContext();
  }
}

/*
  Wraps the function argument to it statements with an async test that runs the
  given work generator context wrapped in act calls.

  The inner generator is called synchronously inside of act.  In order to wait for
  promises that are also wrapped in act, simply yield* runPromise(promise);
 */
export function runInActs(work: () => Generator<any, any, any>) {
  return async function () {
    const context = new Cancellable();
    try {
      let done: boolean | undefined;
      let value: any;

      await context.run(function* () {
        const gen = work();
        act(() => {
          ({done, value} = gen.next())
        });
        while (!done) {
          let resolved: any;
          try {
            if (value instanceof Function) {
              act(() => {
                (resolved = value())
              });
            } else {
              console.log('integration going...')
              yield act(() => Promise.resolve(value).then(v => {
                console.log('resolved', resolved);
                resolved = v
              }));
            }
          } catch (e) {
            act(() => {
              ({done, value} = gen.throw(e))
            })
            continue;
          }
          act(() => {
            ({done, value} = gen.next(resolved))
          });
        }
        return;
      }());
    } catch (e) {
      console.error('Exception in generator', e)
      throw e;
    } finally {
      console.log('ok cancelling from integration')
      context.cancel();
    }
  }
}