import {workflowsResponse} from "../../test_utils/fixtures/workflows-response";
import {integrateElement, setupBefore} from "../../test_utils/integration";
import {useContext} from "react";
import {VulcanContext} from "../vulcan_context";


describe('useWorkflowsLoading', () => {
  const integrated = setupBefore(integrateElement);
  const getWorkflowsMock = setupBefore(() => integrated.value.blockingAsyncMock('getWorkflows'));
  const contextData = setupBefore(() => integrated.value.runHook(() => useContext(VulcanContext)));

  it('works', async () => {
    const {stateRef} = contextData.value;
    expect(stateRef.current.workflows).toEqual([]);
    await getWorkflowsMock.value.respond(() => {
      return Promise.resolve(workflowsResponse)
    });

    expect(stateRef.current.workflows).toEqual(workflowsResponse.workflows);
    expect(getWorkflowsMock.value.pendingCount()).toEqual(0);
  });
});
