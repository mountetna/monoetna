import * as React from "react";
import {integrateElement} from "../../test_utils/integration";
import {asyncFn, AsyncMock, createFakeStorage} from "../../test_utils/mocks";
import {setStatus, setWorkflow} from "../../actions/vulcan";
import {SessionStatusResponse, Workflow} from "../../api_types";
import {
  createStatusFixture, createStepFixture, createStepStatusFixture, createWorkflowFixture
} from "../../test_utils/fixtures";
import {act} from "react-test-renderer";
import {defaultContext} from "../vulcan_context";
import {WorkflowBuilder} from "../../test_utils/workflow_builder";

describe('useSessionSync', () => {
  describe('polling updates', () => {
    let setup: ReturnType<typeof integrateElement>;
    let workflowBuilder: WorkflowBuilder;
    let postInputsMock: AsyncMock<SessionStatusResponse>;
    let pollStatusMock: AsyncMock<SessionStatusResponse>;

    beforeEach(async () => {
      postInputsMock = asyncFn<SessionStatusResponse>()
      pollStatusMock = asyncFn<SessionStatusResponse>()

      setup = integrateElement(() => null, {
        providerOverrides: {
          pollStatus: pollStatusMock.jestMock,
          postInputs: postInputsMock.jestMock,
          logActions: true,
        }
      });

      workflowBuilder = WorkflowBuilder.fromSetup(setup);

      await workflowBuilder.setWorkflow('test-workflow');
      await workflowBuilder.addStep('testStep');
      await workflowBuilder.setStatus('testStep', 'running');
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

    });

    describe('when there is no running step status', () => {
      beforeEach(async () => {
        await workflowBuilder.setStatus('testStep', 'error');
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