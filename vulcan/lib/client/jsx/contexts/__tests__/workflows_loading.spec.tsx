import {act} from "react-test-renderer";
import * as React from "react";
import {delay} from "etna-js/spec/helpers";
import {Workflow, WorkflowsResponse} from "../../api_types";
import {workflowsResponse} from "../../test_utils/fixtures/workflows-response";
import {integrateElement, IntegrateElementParams} from "../../test_utils/integration";
import {asyncFn, AsyncMock} from "../../test_utils/mocks";


describe('useWorkflowsLoading', () => {
  let setup: ReturnType<typeof integrateElement>;
  let getWorkflowsMock: AsyncMock<WorkflowsResponse>;
  let overrides: IntegrateElementParams['providerOverrides'];

  beforeEach(() => {
    getWorkflowsMock = asyncFn();
    overrides = {
      getWorkflows: getWorkflowsMock.jestMock,
    };

    setup = integrateElement(() => null, { providerOverrides: overrides });
  })

  it('works', async () => {
    const {updateMatching, contextData, replaceOverrides} = setup;

    await act(async function () {
      await getWorkflowsMock.awaitCall(false);
      await delay(100); // Ensure there isn't some nasty loop
      expect(getWorkflowsMock.hasPendingRequest()).toBeFalsy()

      getWorkflowsMock.reset();
      // Trigger an empty update.
      replaceOverrides({...overrides});
      await delay(100); // Shouldn't cause any further calls
      expect(getWorkflowsMock.hasPendingRequest()).toBeFalsy()

      // Substantive updates to params does reload.
      getWorkflowsMock.reset();
      replaceOverrides({ ...overrides, params: {a: 2} })
      const [[resolve]] = await getWorkflowsMock.awaitCall(false);
      expect(contextData.state.workflows).not.toBe(workflowsResponse);
      resolve(workflowsResponse);
      await updateMatching(() => contextData.state.workflows === workflowsResponse.workflows);
    });
  });
});
