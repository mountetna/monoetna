import {awaitBefore, integrateElement, setupBefore} from "../../test_utils/integration";
import {useWorkflowUtils} from "../../test_utils/workflow_utils";
import {useContext} from "react";
import {VulcanContext} from "../vulcan_context";
import {addValidationErrors, setBufferedInput, setInputs} from "../../actions/vulcan_actions";

describe('useDataBuffering', () => {
  beforeEach(() => {
    jest.useFakeTimers();
  })

  afterEach(() => {
    jest.useRealTimers();
  })

  const integrated = setupBefore(integrateElement);
  const pollStatusMock = setupBefore(() => integrated.value.blockingAsyncMock('pollStatus'));
  const workflowHelpers = setupBefore(() => integrated.value.runHook(() => useWorkflowUtils()));
  const contextData = setupBefore(() => integrated.value.runHook(() => useContext(VulcanContext)));

  const setupWorkflow = awaitBefore(async () => {
    workflowHelpers.value.setWorkflow('test');
  })

  const setupPrimaryInputs = awaitBefore(async () => {
    workflowHelpers.value.addPrimaryInput('a');
    workflowHelpers.value.addPrimaryInput('b');
    workflowHelpers.value.addPrimaryInput('c');
  })

  const aStep = awaitBefore(async() => {
    return workflowHelpers.value.addStep('astep', {
      out: ['a', 'b'],
      run: 'ui-queries/some-input',
    })
  })

  const bStep = awaitBefore(async() => {
    return workflowHelpers.value.addStep('bstep', {
      out: ['a', 'b'],
      run: 'ui-queries/some-input',
    })
  })

  const setupInputs = awaitBefore(async () => {
    contextData.value.dispatch(setBufferedInput({
      "astep/a": [true],
      "a": [true],
      "bstep/b": [true],
      "bstep/a": [true],
      "c": [true],
    }))
  })

  // Helps verify that other inputs don't get cleared when updating.
  const addOtherInput = awaitBefore(async () => {
    contextData.value.dispatch(setInputs({"v": 2}));
  })

  describe('commitSessionInputs', () => {
    const committedStepName = setupBefore(() => null as string | null);

    const doCommit = awaitBefore(async () => {
      const {commitSessionInputChanges} = contextData.value;
      commitSessionInputChanges(committedStepName.value);
    });

    describe('for primary inputs', () => {
      it('pushes those into the session.inputs, then executes a poll', async () => {
        const {stateRef} = contextData.value;
        expect(Object.keys(stateRef.current.bufferedInputValues)).toEqual(['astep/a', 'bstep/b', 'bstep/a'])
        expect(stateRef.current.session.inputs).toEqual({
          "a": true,
          "c": true,
          "v": 2,
        });
        expect(pollStatusMock.value.pendingCount()).toEqual(1);
        const [[requestedSession]] = await pollStatusMock.value.channel.receive();
        expect(requestedSession.inputs).toEqual(stateRef.current.session.inputs);
      })
    })

    describe('for a step', () => {
      committedStepName.replace(() => 'astep');

      it('pushes those into the session.inputs, then executes a poll', async () => {
        const {stateRef} = contextData.value;
        expect(Object.keys(stateRef.current.bufferedInputValues)).toEqual(['a', 'bstep/b', 'bstep/a', 'c'])
        expect(stateRef.current.session.inputs).toEqual({
          "astep/a": true,
          "v": 2,
        });
        expect(pollStatusMock.value.pendingCount()).toEqual(1);
        const [[requestedSession]] = await pollStatusMock.value.channel.receive();
        expect(requestedSession.inputs).toEqual(stateRef.current.session.inputs);
      })

      describe('with unrelated validation errors', () => {
        const addErrors = awaitBefore(async () => {
          contextData.value.dispatch(addValidationErrors('bstep', 'the label', ['some error']))
        })

        doCommit.replace(async factory => {
          await addErrors.ensure();
          return factory();
        });

        it('still pushes to session inputs and polls', () => {
          const {stateRef} = contextData.value;
          expect(stateRef.current.session.inputs).toEqual({
            "astep/a": true,
            "v": 2,
          });
          expect(pollStatusMock.value.pendingCount()).toEqual(1);
        })

        describe('and related validation errors', () => {
          addErrors.replace(async factory => {
            await factory();
            contextData.value.dispatch(addValidationErrors('astep', 'the label', ['some error 2']))
          })

          it('does not poll or push to session inputs', () => {
            const {stateRef} = contextData.value;
            expect(stateRef.current.session.inputs).toEqual({
              "v": 2,
            });
            expect(pollStatusMock.value.pendingCount()).toEqual(0);
          })
        })
      })
    })
  });

  describe('cancelInputChanges', () => {
    const cancelledStepName = setupBefore(() => null as string | null);

    const doCancel = awaitBefore(async () => {
      const {cancelInputChanges} = contextData.value;
      cancelInputChanges(cancelledStepName.value);
    })

    describe('for primary inputs', () => {
      it('removes only and all outputs from primary inputs', async () => {
        const {stateRef} = contextData.value;
        expect(Object.keys(stateRef.current.bufferedInputValues)).toEqual(['astep/a', 'bstep/b', 'bstep/a'])
      })
    })

    describe('for a step', () => {
      cancelledStepName.replace(() => 'bstep');

      it('removes only and all outputs from that step', async () => {
        const {stateRef} = contextData.value;
        expect(Object.keys(stateRef.current.bufferedInputValues)).toEqual(['astep/a', 'a', 'c'])
      })
    })
  })
});