import {Provider} from "react-redux";
import * as React from "react";
import 'regenerator-runtime';
import {VulcanStore} from "./vulcan_store";
import {defaultContext, ProviderProps, VulcanContextData, VulcanProvider} from "./contexts/vulcan_context";
import {
  defaultSessionStatusResponse,
  defaultWorkflow,
} from "./api_types";
import {workflowsResponse} from "./test_utils/fixtures/workflows-response";
import {defaultSession, defaultVulcanState} from "./reducers/vulcan_reducer";
import {ReactElement} from "react";

type PromiseInner<T extends Promise<any>> = T extends Promise<infer R> ? R : any;
function ioc<F extends (...p: any[]) => Promise<any>, R = PromiseInner<ReturnType<F>>>(
  g: () => Generator<unknown, unknown, [Parameters<F>, (r: R) => void, (e: any) => void]>,
): (...p: Parameters<F>) => Promise<R> {
  const gen = g();
  return (...p: Parameters<F>) => {
    return new Promise<R>((resolve, reject) => {
      gen.next([p, resolve, reject]);
    })
  };
}

function handle<R, P extends any[]>(request: [P, (v: R) => void, (e: any) => void], f: (...p: P) => Promise<R>) {
  const [params, resolve, reject] = request;
  f(...params).then(resolve, reject);
}

export const StorybookMockDecorator = (Story: () => ReactElement) => {
  const store = VulcanStore();

  return <Provider store={store}>
    <VulcanProvider {...mocks}>
      <Story/>
    </VulcanProvider>
  </Provider>;
}


// @ts-ignore
global.ROUTES = {
  workflow: () => `/workflow`,
  workflow_vignette: (workflow_name: string) => `/workflow/${workflow_name}/vignette`,
  fetch_workflows: () => `/api/workflows`,
  submit: (project_name: string, workflow_name: string) =>
    `/api/${project_name}/session/${workflow_name}`,
  status: (project_name: string, workflow_name: string) =>
    `/api/${project_name}/session/${workflow_name}/status`,
};

// @ts-ignore
global.CONFIG = {
  vulcan_host: "",
}

const mocks: Partial<VulcanContextData> & Partial<ProviderProps> = {
  // logActions: true,
  pollStatus: ioc<typeof defaultContext.pollStatus>(function* () {
    while (true) {
      handle(yield, async (session) => {
        return defaultSessionStatusResponse;
      })
    }
  }),
  postInputs: ioc<typeof defaultContext.postInputs>(function*() {
    while (true) {
      handle(yield, async (session) => {
        return defaultSessionStatusResponse;
      })
    }
  }),
  getData: ioc<typeof defaultContext.getData>(function*() {
    while (true) {
      handle(yield, async (url) => {
        return "";
      })
    }
  }),
  getWorkflows: ioc<typeof defaultContext.getWorkflows>(function* () {
  }),
};

