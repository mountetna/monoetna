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
import {createContext, ReactElement, useState} from "react";

type Overrides = Partial<VulcanContextData> & Partial<ProviderProps>;

const defaultOverridesContext = {
  setOverrides(overrides: Overrides) {},
}

export const MockContext = createContext(defaultOverridesContext);

export const StorybookMockDecorator = (Story: () => ReactElement) => {
  const store = VulcanStore();
  const [overrides, setOverrides] = useState({} as Overrides);

  return <Provider store={store}>
    <MockContext.Provider value={{setOverrides}}>
      <VulcanProvider {...{...mocks, ...overrides}}>
        <Story/>
      </VulcanProvider>
    </MockContext.Provider>
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

const {
  getWorkflows,
  getData,
  pollStatus,
  postInputs,
} = defaultContext;

const mocks: Overrides = {
  // logActions: true,
  getWorkflows,
  getData,
  pollStatus,
  postInputs,
};

