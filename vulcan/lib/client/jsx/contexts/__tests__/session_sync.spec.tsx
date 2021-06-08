import * as React from "react";
import {integrateElement} from "../../test_utils/integration";
import {asyncFn, AsyncMock, createFakeStorage} from "../../test_utils/mocks";
import {setStatus, setWorkflow, setWorkflows} from "../../actions/vulcan";
import {defaultSessionStatusResponse, SessionStatusResponse, StatusString, StepStatus, Workflow} from "../../api_types";
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

    async function setupStablePollingState(status: StatusString) {
      const {dispatch, contextData} = setup;
      await act(async () => {
        await dispatch(setStatus(createStatusFixture(workflow, createStepStatusFixture({name: "testStep", status}))));
        const [[resolve]] = await pollStatusMock.awaitCall(false);
        resolve({ ...defaultSessionStatusResponse, status: contextData.stateRef.current.status, session: contextData.stateRef.current.session });
      });
    }

    beforeEach(async () => {
      postInputsMock = asyncFn<SessionStatusResponse>()
      pollStatusMock = asyncFn<SessionStatusResponse>()

      const storage = createFakeStorage();
      workflow = createWorkflowFixture({
        steps: [
          [createStepFixture({ name: "testStep" })],
        ]
      });

      setup = integrateElement(() => null, {
        providerOverrides: {
          storage,
          getWorkflows: defaultContext.getWorkflows,
          pollStatus: pollStatusMock.jestMock,
          postInputs: postInputsMock.jestMock,
        }
      });

      const {dispatch} = setup;
      await act(async () => {
        await dispatch(setWorkflow(workflow, 'test'));
        await setupStablePollingState("running");
      });
    })

    it('polls for a status update', async () => {
      const {contextData} = setup;
      const call = pollStatusMock.awaitCall(true);

      await act(async () => {
        contextData.requestPoll();
      });

      const [_, args] = await call;
      expect(args[0]).toBe(contextData.stateRef.current.session)
    })

    describe('when there is no running step status', () => {
      beforeEach(async () => {
        await setupStablePollingState("complete");
      })

      it('does not continue to poll', async () => {
        const {contextData} = setup;
        const call = pollStatusMock.awaitCall(true);

        await act(async () => {
          contextData.requestPoll();
        });

        expect(pollStatusMock.hasPendingRequest()).toBeFalsy();
        expect(postInputsMock.hasPendingRequest()).toBeFalsy();
      });

      describe('but the request is a post', () => {
        it('does post', async () => {
          const {contextData} = setup;
          const call = pollStatusMock.awaitCall(true);

          await act(async () => {
            contextData.requestPoll(true);
          });

          expect(pollStatusMock.hasPendingRequest()).toBeFalsy();
          expect(postInputsMock.hasPendingRequest()).toBeTruthy();
        })
      })
    })
  })
})