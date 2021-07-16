import * as React from "react";
import {integrateElement, runInActs} from "../../test_utils/integration";
import {countIter, makeBlockingAsyncMock} from "../../test_utils/mocks";
import {
  createStepStatusFixture, createUpdatedStatusFixture,
  createUpdatedStatusResponseFixture,
} from "../../test_utils/fixtures";
import {act} from "react-test-renderer";
import {defaultContext} from "../vulcan_context";
import {runPromise} from "etna-js/utils/cancellable_helpers";
import {useWorkflowUtils} from "../../test_utils/workflow_utils";

describe('useSessionSync', () => {
  beforeEach(() => {
    jest.useFakeTimers();
  })

  afterEach(() => {
    jest.useRealTimers();
  })

  function* runGen<T>(gen: Generator<any, T, any>) {
    return yield* gen;
  }

  function* setupPollingTest() {
    const {update, runHook} = integrateElement();
    const pollStatusMock = makeBlockingAsyncMock(defaultContext.pollStatus);
    const postInputsMock = makeBlockingAsyncMock(defaultContext.postInputs);
    const workflowHelpers = yield* runHook(() => useWorkflowUtils());
    workflowHelpers.setWorkflow('test');

    const {requestPoll, stateRef} = yield* update({ updateOverrides: curOverrides => ({
      ...curOverrides, pollStatus: pollStatusMock.mock, postInputs: postInputsMock.mock
    })});
    return {pollStatusMock, requestPoll, stateRef, workflowHelpers};
  }

  describe('requestPoll', () => {
    it('polls for a status update', runInActs(function* () {
      const {pollStatusMock, requestPoll, stateRef} = yield* runGen(setupPollingTest());

      requestPoll();
      console.log('running promise...')
      const [args, respond] = yield* runPromise(pollStatusMock.channel.receive());
      console.log('ok')
      expect(args[0]).toBe(stateRef.current.session)

      yield jest.advanceTimersByTime(50000);
      expect(countIter(pollStatusMock.channel.pending())).toEqual(0);
    }))

    describe('ongoing polling', () => {
      it('polls every 3 seconds until there remain no running steps', runInActs(function* () {
        const {getContext, replaceOverrides} = integrateElement();
        const pollStatusMock = makeBlockingAsyncMock(defaultContext.pollStatus);
        const postInputsMock = makeBlockingAsyncMock(defaultContext.postInputs);
        const {requestPoll, stateRef} = yield* getContext();

        yield;
        requestPoll(true);

        const [_, respond] = yield* runPromise(postInputsMock.channel.receive());
        respond.resolve(createUpdatedStatusResponseFixture(stateRef.current, {}));

        yield jest.advanceTimersByTime(100);

        expect(stateRef.current.pollingState).toBeTruthy();
        expect(pollStatusMock.hasPendingRequest()).toBeFalsy();

        yield jest.advanceTimersByTime(1000);

        // ([[resolve]] = await pollStatusMock.awaitCall(false));
        //
        // await resolve(createUpdatedStatusResponseFixture(contextData.state, {}));
        // await act(async () => { jest.advanceTimersByTime(1000); });
        //
        // expect(contextData.state.pollingState).toBeTruthy();
        //
        // ([[resolve]] = await pollStatusMock.awaitCall(false));
        // await resolve(createUpdatedStatusResponseFixture(contextData.state, {
        //   status: createUpdatedStatusFixture(workflowUtils.workflow, workflowUtils.status,
        //     createStepStatusFixture({ name: 'testStep', status: 'complete' }))
        // }));
        //
        // expect(pollStatusMock.hasPendingRequest()).toBeFalsy();
        //
        // await act(async () => {
        //   jest.advanceTimersByTime(1000);
        // });
        //
        // expect(pollStatusMock.hasPendingRequest()).toBeFalsy();
        // expect(contextData.state.pollingState).toBeFalsy();
      }))
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