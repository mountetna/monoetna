import React from 'react';
import {awaitBefore, integrateElement, setupBefore} from "../../test_utils/integration";
import {useWorkflowUtils} from "../../test_utils/workflow_utils";
import {useContext} from "react";
import {VulcanContext} from "../vulcan_context";
import {addValidationErrors, setBufferedInput, setInputs} from "../../actions/vulcan_actions";
import {BufferedInputsContext, WithBufferedInputs} from "../input_state_management";
import {PRIMARY_INPUTS} from '../../api_types';

describe('useDataBuffering', () => {
  beforeEach(() => {
    jest.useFakeTimers();
  })

  afterEach(() => {
    jest.useRealTimers();
  })

  const stepName = setupBefore(() => '' as string);
  const integrated = setupBefore(() => integrateElement((hook, {dispatch, commitSessionInputChanges}) =>
    <WithBufferedInputs stepName={stepName.value} dispatch={dispatch}
                        commitSessionInputChanges={commitSessionInputChanges}>
      {hook}
    </WithBufferedInputs>));
  const pollStatusMock = setupBefore(() => integrated.value.blockingAsyncMock('pollStatus'));
  const workflowHelpers = setupBefore(() => integrated.value.runHook(() => useWorkflowUtils()));
  const contextData = setupBefore(() => integrated.value.runHook(() => useContext(VulcanContext)));
  const bufferInputsContextData = setupBefore(() => integrated.value.runHook(() => useContext(BufferedInputsContext)))

  const setupWorkflow = awaitBefore(async () => {
    workflowHelpers.value.setWorkflow('test');
  })

  const setupPrimaryInputs = awaitBefore(async () => {
    workflowHelpers.value.addPrimaryInput('a');
    workflowHelpers.value.addPrimaryInput('b');
    workflowHelpers.value.addPrimaryInput('c');
  })

  const aStep = awaitBefore(async () => {
    return workflowHelpers.value.addStep('astep', {
      out: ['a', 'b'], run: 'ui-queries/some-input',
    })
  })

  const bStep = awaitBefore(async () => {
    return workflowHelpers.value.addStep('bstep', {
      out: ['a', 'b'], run: 'ui-queries/some-input',
    })
  })

  const setupInputs = awaitBefore(async () => {
    bufferInputsContextData.value.setInputs({
      "astep/a": [true], "a": [true], "bstep/b": [true], "bstep/a": [true], "c": [true],
    })
  })

  awaitBefore(async () => {
    expect([...contextData.value.stateRef.current.bufferedSteps].sort()).toEqual([
      stepName.value
    ])
  })

  // Helps verify that other inputs don't get cleared when updating.
  const addOtherInput = awaitBefore(async () => {
    contextData.value.dispatch(setInputs({"v": 2}));
  })

  describe('commitSessionInputs', () => {
    const doCommit = awaitBefore(async () => {
      bufferInputsContextData.value.commitInputs();
    });

    describe('for primary inputs', () => {
      stepName.replace(() => PRIMARY_INPUTS);

      it('pushes those into the session.inputs, then executes a poll', async () => {
        const {stateRef} = contextData.value;
        const {inputs} = bufferInputsContextData.value;
        expect(Object.keys(inputs)).toEqual([])
        expect(stateRef.current.validationErrors).toEqual([]);
        expect(stateRef.current.session.inputs).toEqual({
          "a": true, "c": true, "v": 2,
        });
        expect(pollStatusMock.value.pendingCount()).toEqual(1);
        const [[requestedSession]] = await pollStatusMock.value.channel.receive();
        expect(requestedSession.inputs).toEqual(stateRef.current.session.inputs);
      })
    })

    describe('for a step', () => {
      stepName.replace(() => 'astep');

      it('pushes those into the session.inputs, then executes a poll', async () => {
        const {stateRef} = contextData.value;
        const {inputs} = bufferInputsContextData.value;
        expect(Object.keys(inputs)).toEqual([])
        expect(stateRef.current.session.inputs).toEqual({
          "astep/a": true, "v": 2,
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
            "astep/a": true, "v": 2,
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
    const doCancel = awaitBefore(async () => {
      bufferInputsContextData.value.cancelInputs();
    })

    it('cleans the buffered inputs and does not effect the state inputs', async () => {
      const {stateRef} = contextData.value;
      const {inputs} = bufferInputsContextData.value;
      expect(Object.keys(inputs)).toEqual([])
      expect(Object.keys(stateRef.current.session.inputs)).toEqual(['v'])
      expect(stateRef.current.bufferedSteps).toEqual([])
    })
  })
});