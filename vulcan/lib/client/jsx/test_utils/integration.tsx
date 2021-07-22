import {ReactElement, useContext} from "react";
import {act, create} from "react-test-renderer";
import {
  defaultContext, ProviderProps, VulcanContext, VulcanContextData, VulcanProvider
} from "../contexts/vulcan_context";
import {Provider} from "react-redux";
import {VulcanStore} from "../vulcan_store";
import * as React from "react";
import {BlockingAsyncMock, createFakeStorage, UnbufferedAsyncMock} from "./mocks";
import {Maybe, some} from "../selectors/maybe";

type InnerElementConstructor = (hookElement: ReactElement, contextData: VulcanContextData) => ReactElement;

export function integrateElement(
  element: ReactElement | null | InnerElementConstructor = null,
  providerOverrides: Partial<ProviderProps & VulcanContextData> = {},
) {
  // Default provide mocks should proceed the spread below.
  providerOverrides = {
    getWorkflows: defaultContext.getWorkflows,
    getData: defaultContext.getData,
    pollStatus: defaultContext.pollStatus,
    postInputs: defaultContext.postInputs,
    storage: createFakeStorage(), ...providerOverrides,
  }

  const store = VulcanStore();

  let preHook: Function | undefined = undefined;
  let lastPreHook: Function | undefined = preHook;
  let prehookResult = [undefined] as [any];
  let renderIdx = 0;
  let PreHookContainer = function PreHookContainer() {
    if (preHook) prehookResult[0] = preHook();
    return null;
  }

  const node = create(rerender())

  return {node, runHook, provideOverrides, replaceElement, blockingAsyncMock, unbufferedAsyncMock};

  // Regenerates the react element tree to create or update.
  // Manages & injects the prehook as well as the inner element.
  function rerender() {
    if (preHook !== lastPreHook) {
      renderIdx++;
      lastPreHook = preHook;
      prehookResult[0] = undefined;
      PreHookContainer = function PreHookContainer() {
        if (preHook) prehookResult[0] = preHook();
        return null;
      }
    }

    const prehook = <PreHookContainer key={renderIdx + "-prehook"}/>;
    const inner = element instanceof Function ?
      <ProvidesElementFunction element={element} prehook={prehook}/> :
      element;

    return <Provider store={store}>
      <VulcanProvider {...providerOverrides}>
        {inner}
      </VulcanProvider>
    </Provider>;
  }

  function ProvidesElementFunction({element, prehook}: { element: InnerElementConstructor, prehook: ReactElement }) {
    const contextData = useContext(VulcanContext);
    return element(prehook, contextData);
  }

  function runHook<T>(hook: () => T): T {
    act(function () {
      preHook = hook;
      node.update(rerender());
    });

    return prehookResult[0] as T;
  }

  function getContext(): VulcanContextData {
    return runHook(function () {
      return useContext(VulcanContext);
    })
  }

  function provideOverrides(p: Partial<ProviderProps & VulcanContextData>) {
    Object.assign(providerOverrides, p);
    return getContext();
  }

  function replaceElement(e: ReactElement | null) {
    element = e;
    return getContext();
  }

  function blockingAsyncMock<K extends keyof VulcanContextData, T>(k: K): VulcanContextData[K] extends (...a: any[]) => Promise<T> ? BlockingAsyncMock<any[], T> : never {
    const mock = new BlockingAsyncMock<any[], T>(defaultContext[k] as any, act);
    provideOverrides({[k]: mock.mock});
    return mock as any;
  }

  function unbufferedAsyncMock<K extends keyof VulcanContextData, T>(k: K): VulcanContextData[K] extends (...a: any[]) => Promise<T> ? UnbufferedAsyncMock<any[], T> : never {
    const mock = new UnbufferedAsyncMock<any[], T>(defaultContext[k] as any, act);
    provideOverrides({[k]: mock.mock});
    return mock as any;
  }
}

export class ValueCell<ValueType> {
  cached: Maybe<ValueType> = null;
  originalFactory: () => Promise<ValueType> | ValueType;

  constructor(
    protected factory: () => Promise<ValueType> | ValueType,
    protected scheduler: (f: () => Promise<void>) => Promise<void> = async f => await f(),
  ) {
    this.originalFactory = factory;
    this.setup();
  }

  get value(): ValueType {
    if (this.cached) return this.cached[0];
    throw new Error('.value not ready, did you forget to await .ensure() in a definition?')
  }

  protected setup() {
    beforeAll(() => {
      this.factory = this.originalFactory;
    })

    beforeEach(() => this.ensure());

    afterEach(() => {
      this.cached = null;
    })
  }

  async ensure() {
    if (!this.cached) {
      await this.scheduler(async () => {
        const v = this.factory();
        this.cached = some(await v);
      });
    }

    if (this.cached) return this.cached[0];
    throw new Error('Ensure did not ensure value creation, bug.');
  }

  replace(newF: (f: () => Promise<ValueType> | ValueType) => (Promise<ValueType> | ValueType)) {
    let lastFactory: () => Promise<ValueType> | ValueType;
    beforeAll(() => {
      const curFactory = lastFactory = this.factory;
      this.factory = () => newF(curFactory)
    })

    afterAll(() => {
      this.factory = lastFactory;
    })
  }
}

export function setupBefore<T>(f: () => T) {
  return new ValueCell(f, async inner => {
    // Run the inner, synchronous
    await inner();
    await act(async () => {
    });
  });
}

export function awaitBefore<T>(f: () => Promise<T>) {
  return new ValueCell(f, act);
}