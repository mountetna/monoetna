import React from 'react';
import {
  awaitBefore,
  integrateElement,
  setupBefore
} from '../../test_utils/integration';
import {useWorkflowUtils} from '../../test_utils/workflow_utils';
import {useContext} from 'react';
import {VulcanContext} from '../vulcan_context';
import {
  addValidationErrors,
  setAutoPassStep,
  setSession,
} from '../../actions/vulcan_actions';
import {
  BufferedInputsContext,
  WithBufferedInputs
} from '../input_state_management';
import {createSessionFixture, createStatusFixture, createStepStatusFixture} from '../../test_utils/fixtures';
import { Switch } from '@material-ui/core';
import { act } from 'react-test-renderer';

describe('useDataBuffering', () => {
  beforeEach(() => {
    jest.useFakeTimers();
  });

  afterEach(() => {
    jest.useRealTimers();
  });

  const stepName = setupBefore(() => 'bstep' as string | null);
  const integrated = setupBefore(() =>
    integrateElement((hook, {dispatch, commitSessionInputChanges}) => (
      <WithBufferedInputs
        stepName={stepName.value}
        dispatch={dispatch}
        commitSessionInputChanges={commitSessionInputChanges}
      >
        {hook}
      </WithBufferedInputs>
    ))
  );
  const pollStatusMock = setupBefore(() =>
    integrated.value.blockingAsyncMock('pollStatus')
  );
  const workflowHelpers = setupBefore(() =>
    integrated.value.runHook(() => useWorkflowUtils())
  );
  const contextData = setupBefore(() =>
    integrated.value.runHook(() => useContext(VulcanContext))
  );
  const bufferInputsContextData = setupBefore(() =>
    integrated.value.runHook(() => useContext(BufferedInputsContext))
  );

  const setupWorkflow = awaitBefore(async () => {
    workflowHelpers.value.setWorkflow('test');
  });

  const testSession = setupBefore(() =>
    createSessionFixture('test', {project_name: 'test'})
  );

  const setupSession = awaitBefore(async () => {
    await testSession.ensure();
    contextData.value.dispatch(setSession(testSession.value));
  });

  const setupPrimaryInputs = awaitBefore(async () => {
    workflowHelpers.value.addPrimaryInput('a');
    workflowHelpers.value.addPrimaryInput('b');
    workflowHelpers.value.addPrimaryInput('c');
  });

  const aStep = awaitBefore(async () => {
    return workflowHelpers.value.addStep('astep', {
      out: ['a', 'b'],
      doc: "blah blah blah",
      run: 'ui-queries/some-input'
    });
  });

  const bStep = awaitBefore(async () => {
    return workflowHelpers.value.addStep('bstep', {
      out: ['a', 'b'],
      doc: "SKIPPABLE. blah blah blah",
      run: 'ui-queries/some-input'
    });
  });

  describe('autoPass', () => {

    describe('for primary inputs', () => {
      stepName.replace(() => null)
      describe('in workflow w/ vignette trigger', () => {
        const addVignetteWithTrigger = awaitBefore(async () => {
          workflowHelpers.value.setVignette('blah blah Primary inputs are skippable blah blah');
        });
        describe('no session inputs yet', () => {
          const setPrimaryInputs = awaitBefore(async () => {
            bufferInputsContextData.value.setInputs(() => {
              return { a: [true], b: [true], c: [true] }
            })
          })
          it('autoPass initiated', async () => {
            const {stateRef} = contextData.value;
            expect(stateRef.current.workflow?.vignette?.includes("Primary inputs are skippable")).toBeTruthy()
            expect(stateRef.current.triggerRun).toEqual([null])
          })
        })
        // describe('with previous session inputs', () => {
        //   testSession.replace(() =>
        //     createSessionFixture('test', {project_name: 'test', inputs: {a: [true]}})
        //   )
        //   const setupSession = awaitBefore(async () => {
        //     await testSession.ensure();
        //     contextData.value.dispatch(setSession(testSession.value));
        //   });
        //   const bufferInputsContextData = setupBefore(() =>
        //     integrated.value.runHook(() => useContext(BufferedInputsContext))
        //   );
        //   describe('before primaries set', () => {
        //     const setPrimaryInputs = awaitBefore(async () => {
        //       bufferInputsContextData.value.setInputs(() => {
        //         return { a: [true], b: [true], c: [true] }
        //       })
        //     })
        //     it('does NOT initiate autoPass', () => {
        //       const {state} = contextData.value;
        //       console.log(state.session.inputs)
        //       console.log("after?")
        //       const {stateRef} = contextData.value;
        //       expect(stateRef.current.workflow?.vignette?.includes("Primary inputs are skippable")).toBeTruthy()
        //       expect(stateRef.current.triggerRun).toEqual([])
        //     })
        //   })
        // })
      })
      describe('for primary inputs, w/out vignette string', () => {
        const addVignetteWithoutTrigger = awaitBefore(async () => {
          workflowHelpers.value.setVignette('blah blah blah blah');
        });
        const setPrimaryInputs = awaitBefore(async () => {
          bufferInputsContextData.value.setInputs(() => {
            return { a: [true], b: [true], c: [true] }
          })
        })
        it('TryCommitThenRun not initiated', () => {
          const {stateRef} = contextData.value;
          expect(stateRef.current.workflow?.vignette?.includes("Primary inputs are skippable")).toBeFalsy()
          expect(stateRef.current.triggerRun).toEqual([])
        })
      });
    })

    describe('for stepUserInputs w/ doc string', () => {
      stepName.replace(() => 'bstep')
      async function pushStepUI() {
        await act(async () => {
          bufferInputsContextData.value.setInputs(() => {
            return { 'bstep/b': [true], 'bstep/a': [true] }
          })
        });
      }
      describe('eligible', () => {
        const setStepInputs = awaitBefore(async () => {
          bufferInputsContextData.value.setInputs(() => {
            return { 'bstep/b': [true], 'bstep/a': [true] }
          })
        })
        it('switch shown', () => {
          expect(contextData.value.stateRef.current.autoPassSteps).toEqual([])
          expect(contextData.value.stateRef.current.tryCommitThenRun).toEqual([])
          expect(contextData.value.stateRef.current.triggerRun).toEqual([])
          const {node} = integrated.value;
          expect(node.root.findAllByType(Switch).length).toEqual(1);
          expect(node.root.findAllByProps({className:'reset-or-commit-inputs'}).length).toEqual(1)
        })
      })
      // fit('adds to autoPassSteps on toggle', async () => {
      //   const {node} = integrated.value;
      //   await act(async () => {
      //     node.root.findByType(Switch).props.onChange();
      //   });
      //   // Why doesn't this update?
      //   console.log(contextData.value.stateRef.current)
      //   expect(contextData.value.stateRef.current.autoPassSteps).toEqual(['bstep'])
      // })
      describe('turned on', () => {
        const setAutoPass = awaitBefore(async () => {
          contextData.value.dispatch(setAutoPassStep(stepName.value));
        })
        describe('no validation errors', () => {
          it('gets through TryCommitThenRun & triggers Run', async () => {
            const {stateRef} = contextData.value;
            const {inputs} = bufferInputsContextData.value;
            expect(stateRef.current.validationErrors).toEqual([]);
            expect(stateRef.current.autoPassSteps).toEqual(['bstep']);
            expect(stateRef.current.tryCommitThenRun).toEqual([]);
            expect(stateRef.current.triggerRun).toEqual([]);
            // Triggered when setting inputs, with tryCommitThenRun also acted on within this await
            await pushStepUI()
            expect(stateRef.current.autoPassSteps).toEqual(['bstep']);
            expect(stateRef.current.tryCommitThenRun).toEqual([]);
            expect(stateRef.current.triggerRun).toEqual(['bstep']);
          })
          // describe('with previous session inputs for step', () => {
          //   testSession.replace(() =>
          //     createSessionFixture('test', {project_name: 'test', inputs: {bstep/a: [true]}})
          //   )
          //   const setupSession = awaitBefore(async () => {
          //     await testSession.ensure();
          //     contextData.value.dispatch(setSession(testSession.value));
          //   });
          //   const bufferInputsContextData = setupBefore(() =>
          //     integrated.value.runHook(() => useContext(BufferedInputsContext))
          //   );
          //   describe('on passing new inputs', () => {
          //     const setStepInputs = awaitBefore(async () => {
          //       bufferInputsContextData.value.setInputs(() => {
          //         return { 'bstep/b': [true], 'bstep/a': [true] }
          //       })
          //     })
          //     it('does NOT initiate autoPass', () => {
          //       const {stateRef} = contextData.value;
          //       expect(stateRef.current.autoPassSteps).toEqual(['bstep']);
          //       expect(stateRef.current.tryCommitThenRun).toEqual([]);
          //       expect(stateRef.current.triggerRun).toEqual([]);
          //     })
          //   })
          // })
        })
        describe('with validation errors', () => {
          const addValidationError = awaitBefore(async () => {
            contextData.value.dispatch(
              addValidationErrors('bstep', 'the label', ['some error'])
            );
          })
          it('NO triggerRun, ', async () => {
            const {stateRef} = contextData.value;
            const {inputs} = bufferInputsContextData.value;
            expect(stateRef.current.validationErrors).toEqual([["bstep","the label",["some error"]]]);
            // Triggered when setting inputs, with tryCommitThenRun also acted on within this await
            await pushStepUI()
            expect(stateRef.current.autoPassSteps).toEqual(['bstep']);
            expect(stateRef.current.tryCommitThenRun).toEqual([]);
            expect(stateRef.current.triggerRun).toEqual([]);
          })
        })
      })
    })
    describe('for stepUserInputs w/out doc string', () => {
      stepName.replace(() => 'astep')
      describe('not eligible', () => {
        it('switch not shown', async () => {
          await act(async () => {
            bufferInputsContextData.value.setInputs(() => {
              return { 'astep/b': [true], 'astep/a': [true] }
            })
          });
          const {node} = integrated.value;
          expect(node.root.findAllByProps({className:'reset-or-commit-inputs'}).length).toEqual(1)
          expect(node.root.findAllByType(Switch).length).toEqual(0);
        })
      })
    })
  })
});
