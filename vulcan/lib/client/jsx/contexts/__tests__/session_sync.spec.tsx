import * as React from "react";
import {integrateElement} from "../../test_utils/integration";
import {asyncFn, AsyncMock, createFakeStorage} from "../../test_utils/mocks";
import {setStatus, setWorkflow} from "../../actions/vulcan";
import {SessionStatusResponse, Workflow} from "../../api_types";
import {
  createStatusFixture,
  createStepFixture,
  createStepStatusFixture, createUpdatedStatusFixture,
  createUpdatedStatusResponseFixture,
  createWorkflowFixture
} from "../../test_utils/fixtures";
import {act} from "react-test-renderer";
import {defaultContext} from "../vulcan_context";
import {WorkflowUtils} from "../../test_utils/workflow_utils";

describe('useSessionSync', () => {
  describe('polling updates', () => {
    let setup: ReturnType<typeof integrateElement>;
    let workflowUtils: WorkflowUtils;
    let postInputsMock: AsyncMock<SessionStatusResponse>;
    let pollStatusMock: AsyncMock<SessionStatusResponse>;

    beforeEach(async () => {
      postInputsMock = asyncFn<SessionStatusResponse>()
      pollStatusMock = asyncFn<SessionStatusResponse>()

      setup = integrateElement(() => null, {
        providerOverrides: {
          pollStatus: pollStatusMock.jestMock,
          postInputs: postInputsMock.jestMock,
        }
      });

      workflowUtils = WorkflowUtils.fromSetup(setup);

      await workflowUtils.setWorkflow('test-workflow');
      await workflowUtils.addStep('testStep');
      await workflowUtils.setStatus('testStep', 'running');
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

    describe('ongoing polling', () => {
      beforeEach(() => {
        jest.useFakeTimers();
      })

      afterEach(() => {
        jest.useRealTimers();
      })

      it('polls every 10 seconds until there remain no running steps', async () => {
        const {contextData} = setup;

        await act(async () => {
          contextData.requestPoll(true);
        });

        let [[resolve]] = await postInputsMock.awaitCall(false);
        await resolve(createUpdatedStatusResponseFixture(contextData.state, {}));

        await act(async () => {
          jest.advanceTimersByTime(100);
        });

        expect(contextData.state.pollingState).toBeTruthy();
        expect(pollStatusMock.hasPendingRequest()).toBeFalsy();

        await act(async () => {
          jest.advanceTimersByTime(1000);
        });

        ([[resolve]] = await pollStatusMock.awaitCall(false));

        await resolve(createUpdatedStatusResponseFixture(contextData.state, {}));
        await act(async () => { jest.advanceTimersByTime(1000); });

        expect(contextData.state.pollingState).toBeTruthy();

        ([[resolve]] = await pollStatusMock.awaitCall(false));
        await resolve(createUpdatedStatusResponseFixture(contextData.state, {
          status: createUpdatedStatusFixture(workflowUtils.workflow, workflowUtils.status,
            createStepStatusFixture({ name: 'testStep', status: 'complete' }))
        }));

        expect(pollStatusMock.hasPendingRequest()).toBeFalsy();

        await act(async () => {
          jest.advanceTimersByTime(1000);
        });

        expect(pollStatusMock.hasPendingRequest()).toBeFalsy();
        expect(contextData.state.pollingState).toBeFalsy();
      })
    });

    describe('when there is no running step status', () => {
      beforeEach(async () => {
        await workflowUtils.setStatus('testStep', 'error');
      })

      it('does not continue to poll', async () => {
        const {contextData} = setup;

        await act(async () => {
          contextData.requestPoll();
        });

        expect(pollStatusMock.hasPendingRequest()).toBeFalsy();
        expect(postInputsMock.hasPendingRequest()).toBeFalsy();
      });

      describe('but the request is a post', () => {
        it('does post', async () => {
          const {contextData} = setup;

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