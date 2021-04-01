import {integrateElement} from "../../test_utils/integration";
import {SessionStatusResponse, VulcanSession, WorkflowsResponse} from "../../api_types";
import {workflowsResponse} from "../../test_utils/fixtures/workflows-response";
import {act} from "react-test-renderer";
import {setInputs, setWorkflow} from "../../actions/vulcan";
import {createStepFixture, createWorkflowFixture, findWorkflowFromResponse} from "../../test_utils/fixtures";

describe('clear_obsolete_inputs', () => {
  it('works', async () => {
    const {contextData, setInput} = integrateElement(() => null, {
      // contextOverrides: {
      //   getWorkflows(): Promise<WorkflowsResponse> {
      //     return Promise.resolve(workflowsResponse);
      //   },
      //   pollStatus(session: VulcanSession): Promise<SessionStatusResponse> {
      //     return new Promise
      //   }
      // }
    });

    await act(async function() {
      contextData.dispatch(setWorkflow(findWorkflowFromResponse(workflowsResponse, "test_concurrent_workflow.cwl")));
      setInput('a', 10);
      // expect(contextData.state.inputs)
    });
  })
})