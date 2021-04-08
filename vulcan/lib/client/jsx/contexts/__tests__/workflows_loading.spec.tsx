import {act} from "react-test-renderer";
import * as React from "react";
import {delay} from "etna-js/spec/helpers";
import {WorkflowsResponse} from "../../api_types";
import {workflowsResponse} from "../../test_utils/fixtures/workflows-response";
import {integrateElement} from "../../test_utils/integration";


describe('useWorkflowsLoading', () => {
  it('works', async () => {
    let numCalls = 0;
    const overrides = {
        getWorkflows(): Promise<WorkflowsResponse> {
          numCalls += 1;
          return Promise.resolve(workflowsResponse);
        }
    };

    const {contextData, updateMatching, replaceOverrides} = integrateElement(() => null, { providerOverrides: overrides });

    await act(async function () {
      await updateMatching(() => contextData.state.workflows === workflowsResponse.workflows);
      expect(numCalls).toEqual(1);
      await delay(100); // Ensure there isn't some nasty loop
      expect(numCalls).toEqual(1);

      replaceOverrides({});
      await delay(100); // Shouldn't cause any further calls
      expect(numCalls).toEqual(1);

      // Substantive updates to params does reload.
      replaceOverrides({ ...overrides, params: {a: 2} })
      await updateMatching(() => numCalls == 2);
      await updateMatching(() => contextData.state.workflows === workflowsResponse.workflows);
    });
  });
});
