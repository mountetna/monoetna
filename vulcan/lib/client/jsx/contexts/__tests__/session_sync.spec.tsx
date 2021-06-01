import * as React from "react";
import {integrateElement} from "../../test_utils/integration";
import {asyncFn, AsyncMock, createFakeStorage} from "../../test_utils/mocks";
import {setStatus, setWorkflow, setWorkflows} from "../../actions/vulcan";
import {SessionStatusResponse, Workflow} from "../../api_types";
import {
  createStatusFixture,
  createStepFixture, createStepStatusFixture,
  createWorkflowFixture
} from "../../test_utils/fixtures";
import {act} from "react-test-renderer";
import {defaultContext} from "../vulcan_context";

describe('useSessionSync', () => {
  describe('polling updates', () => {
    let setup: ReturnType<typeof integrateElement>;
    let workflow: Workflow;
    let postInputsMock: AsyncMock<SessionStatusResponse>;
    let pollStatusMock: AsyncMock<SessionStatusResponse>;

    beforeEach(async () => {
      postInputsMock = asyncFn<SessionStatusResponse>()
      pollStatusMock = asyncFn<SessionStatusResponse>()

      const storage = createFakeStorage();
      workflow = createWorkflowFixture({
        steps: [
          [createStepFixture({ name: "a" })],
        ]
      });

      setup = integrateElement(() => null, {
        providerOverrides: {
          storage,
          getWorkflows: defaultContext.getWorkflows,
          pollStatus: pollStatusMock[0],
          postInputs: postInputsMock[0],
        }
      });

      const {dispatch} = setup;
      await act(async () => {
        await dispatch(setWorkflow(workflow, 'test'));
        await dispatch(setStatus(createStatusFixture(workflow, createStepStatusFixture({name: "a", status: "running"}))));
      })
    })

    it('polls for a status update', async () => {
      const {contextData} = setup;
      expect(pollStatusMock[0]).toHaveBeenCalledWith(contextData.stateRef.current.session)
      pollStatusMock[0].mockClear();

      await act(async () => {
        contextData.requestPoll();
      });

      expect(pollStatusMock[0]).toHaveBeenCalledWith(contextData.stateRef.current.session)
    })

    describe('when there is no running step status', () => {
      beforeEach(async () => {
        await act(async () => {
          const {dispatch} = setup;
          await dispatch(setStatus(createStatusFixture(workflow, createStepStatusFixture({name: "a", status: "complete"}))));
        })
      })

      it('does not continue to poll', async () => {
        const {contextData} = setup;
        pollStatusMock[0].mockClear();

        await act(async () => {
          contextData.requestPoll();
        });

        expect(pollStatusMock[0]).not.toHaveBeenCalled();
      })

      describe('but the request is a post', () => {
        it('does post', async () => {
          const {contextData} = setup;
          pollStatusMock[0].mockClear();

          await act(async () => {
            contextData.requestPoll(true);
          });

          expect(pollStatusMock[0]).not.toHaveBeenCalled();
          expect(postInputsMock[0]).toHaveBeenCalled();
        })
      })
    })
  })
})