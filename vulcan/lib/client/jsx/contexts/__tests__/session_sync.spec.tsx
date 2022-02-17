import * as React from 'react';
import {
  awaitBefore,
  integrateElement,
  setupBefore
} from '../../test_utils/integration';
import {
  createStepStatusFixture,
  createUpdatedStatusFixture,
  createUpdatedStatusResponseFixture
} from '../../test_utils/fixtures';
import {useWorkflowUtils} from '../../test_utils/workflow_utils';
import {useContext} from 'react';
import {VulcanContext} from '../vulcan_context';
import {act} from 'react-test-renderer';

describe('useSessionSync', () => {
  beforeEach(() => {
    jest.useFakeTimers();
  });

  afterEach(() => {
    jest.useRealTimers();
  });

  const integrated = setupBefore(integrateElement);
  const pollMock = setupBefore(() =>
    integrated.value.blockingAsyncMock('pollStatus')
  );
  const postMock = setupBefore(() =>
    integrated.value.blockingAsyncMock('postInputs')
  );
  const workflowHelpers = setupBefore(() =>
    integrated.value.runHook(() => useWorkflowUtils())
  );
  const contextData = setupBefore(() =>
    integrated.value.runHook(() => useContext(VulcanContext))
  );

  awaitBefore(async () => {
    workflowHelpers.value.setWorkflows([
      {
        name: 'test'
      }
    ]);
    workflowHelpers.value.setWorkflow('test');
    workflowHelpers.value.addStep('astep');
    workflowHelpers.value.setStatus('astep', 'complete');
  });

  const pollStarted = awaitBefore(async () => {
    contextData.value.requestPoll();
  });

  describe('requestPoll', () => {
    it('polls for a status update', async () => {
      const {stateRef} = contextData.value;
      const [args] = await pollMock.value.channel.receive();
      expect(args[0]).toBe(stateRef.current.session);
      jest.advanceTimersByTime(50000);

      expect(pollMock.value.pendingCount()).toEqual(0);
    });

    describe('ongoing polling', () => {
      pollStarted.replace(async () => {
        contextData.value.requestPoll(true);
      });

      it('polls every 3 seconds until there remain no running steps', async () => {
        const {stateRef} = contextData.value;
        console.log(stateRef);
        expect(stateRef.current.pollingState).toBeGreaterThan(0);
        await postMock.value.respond(() =>
          createUpdatedStatusResponseFixture(stateRef.current, {
            status: createUpdatedStatusFixture(
              workflowHelpers.value.workflow,
              workflowHelpers.value.status,
              createStepStatusFixture({name: 'astep', status: 'running'})
            )
          })
        );

        act(() => {
          jest.advanceTimersByTime(100);
        });

        expect(stateRef.current.pollingState).toBeTruthy();
        expect(pollMock.value.pendingCount()).toEqual(0);
        act(() => {
          jest.advanceTimersByTime(3000);
        });

        await pollMock.value.respond(() =>
          createUpdatedStatusResponseFixture(stateRef.current, {})
        );

        act(() => {
          jest.advanceTimersByTime(3000);
        });
        expect(stateRef.current.pollingState).toBeGreaterThan(0);

        await pollMock.value.respond(() =>
          createUpdatedStatusResponseFixture(stateRef.current, {
            status: createUpdatedStatusFixture(
              workflowHelpers.value.workflow,
              workflowHelpers.value.status,
              createStepStatusFixture({name: 'astep', status: 'complete'})
            )
          })
        );

        expect(pollMock.value.pendingCount()).toEqual(0);
        expect(stateRef.current.pollingState).toEqual(0);

        act(() => {
          jest.advanceTimersByTime(5000);
        });
        expect(pollMock.value.pendingCount()).toEqual(0);
      });
    });
  });
});
