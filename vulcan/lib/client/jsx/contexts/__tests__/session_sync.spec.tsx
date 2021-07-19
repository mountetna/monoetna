import * as React from "react";
import {integrateElement, runInActs} from "../../test_utils/integration";
import {countIter, makeBlockingAsyncMock} from "../../test_utils/mocks";
import {
  createStepStatusFixture, createUpdatedStatusFixture,
  createUpdatedStatusResponseFixture,
} from "../../test_utils/fixtures";
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
    workflowHelpers.addStep('astep');
    workflowHelpers.setStatus('astep', 'complete');

    const {requestPoll, stateRef} = yield* update({ updateOverrides: curOverrides => ({
      ...curOverrides, pollStatus: pollStatusMock.mock, postInputs: postInputsMock.mock
    })});
    return {pollStatusMock, postInputsMock, requestPoll, stateRef, workflowHelpers};
  }

  describe('requestPoll', () => {
    it('polls for a status update', runInActs(function* () {
      const {pollStatusMock, requestPoll, stateRef} = yield* runGen(setupPollingTest());

      requestPoll();
      const [args] = yield* runPromise(pollStatusMock.channel.receive());
      expect(args[0]).toBe(stateRef.current.session)

      yield jest.advanceTimersByTime(50000);
      expect(countIter(pollStatusMock.channel.pending())).toEqual(0);
    }))

    describe('ongoing polling', () => {
      it('polls every 3 seconds until there remain no running steps', runInActs(function* () {
        const {workflowHelpers, pollStatusMock, postInputsMock, requestPoll, stateRef} = yield* runGen(setupPollingTest());
        requestPoll(true);

        yield
        expect(stateRef.current.pollingState).toBeGreaterThan(0);
        let [_, respond] = yield* runPromise(postInputsMock.channel.receive());
        yield respond.resolve(createUpdatedStatusResponseFixture(stateRef.current, {
          status: createUpdatedStatusFixture(workflowHelpers.workflow, workflowHelpers.status, createStepStatusFixture(
            {name: 'astep', status: 'running'}
          ))
        }));

        yield jest.advanceTimersByTime(100);

        expect(stateRef.current.pollingState).toBeTruthy();
        expect(countIter(pollStatusMock.channel.pending())).toEqual(0);
        yield jest.advanceTimersByTime(3000);

        ([_, respond] = yield* runPromise(pollStatusMock.channel.receive()));
        yield respond.resolve(createUpdatedStatusResponseFixture(stateRef.current, {}));

        yield jest.advanceTimersByTime(3000);
        ([_, respond] = yield* runPromise(pollStatusMock.channel.receive()));
        expect(stateRef.current.pollingState).toBeGreaterThan(0);
        yield respond.resolve(createUpdatedStatusResponseFixture(stateRef.current, {
          status: createUpdatedStatusFixture(workflowHelpers.workflow, workflowHelpers.status, createStepStatusFixture(
            {name: 'astep', status: 'complete'}
          ))
        }));

        expect(countIter(pollStatusMock.channel.pending())).toEqual(0);
        expect(stateRef.current.pollingState).toEqual(0);

        yield jest.advanceTimersByTime(5000);
        expect(countIter(pollStatusMock.channel.pending())).toEqual(0);
      }))
    });
  })
})