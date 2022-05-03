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
  setBufferedInput,
  setInputs,
  setSession,
  setStatus
} from '../../actions/vulcan_actions';
import {
  BufferedInputsContext,
  WithBufferedInputs
} from '../input_state_management';
import {createSessionFixture, createStatusFixture, createStepStatusFixture} from '../../test_utils/fixtures';
import { Switch } from '@material-ui/core';
import { text } from '../../test_utils/rendered';
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

  const setVignette = awaitBefore(async () => {
    workflowHelpers.value.setVignette('blah blah Primary inputs are skippable blah blah');
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

  const setupInputs = awaitBefore(async () => {
    bufferInputsContextData.value.setInputs({
      'astep/a': [true],
      a: [true],
      b: [true],
      'bstep/b': [true],
      'bstep/a': [true],
      c: [true]
    });
  });

  awaitBefore(async () => {
    expect(
      [...contextData.value.stateRef.current.bufferedSteps].sort()
    ).toEqual([stepName.value]);
  });

  describe('autoPass', () => {
    // const doCommit = awaitBefore(async () => {
    //   bufferInputsContextData.value.commitInputs();
    // });

    describe('for primary inputs', () => {
      const setVignette = awaitBefore(async () => {
        workflowHelpers.value.setVignette('blah blah Primary inputs are skippable blah blah');
      }); 
      //   integrated.replace(() => integrateElement( () => <PrimaryInputs/> ))
      //   fit('TryCommitThenRun initiated', () => {
      //     expect(contextData.value.stateRef.current.workflow?.vignette?.includes("Primary inputs are skippable")).toBeTruthy()
      //     console.log(bufferInputsContextData.value.inputs)
      //     console.log(contextData.value)
      //   })
      //   it('clears TryCommitThenRun & triggers Run if commit is successful', () => {})
      //   it('clears TryCommitThenRun, if commit has validation errors', () => {})
      // })
      // describe('for primary inputs, w/out vignette string', () => {
      //   setVignette.replace(async () => {
      //     workflowHelpers.value.setVignette('blah blah blah blah');
      //   })
      //   it('TryCommitThenRun not initiated', () => {
      //     console.log(bufferInputsContextData.value.inputs)
      //     console.log(contextData.value.stateRef.current.workflow)
      //   })
      // });
    })
    describe('for stepUserInputs w/ doc string', () => {
      
      stepName.replace(() => 'bstep')
      
      it('eligible - switch shown', () => {
        expect(contextData.value.stateRef.current.autoPassSteps).toEqual([])
        expect(contextData.value.stateRef.current.tryCommitThenRun).toEqual([])
        expect(contextData.value.stateRef.current.triggerRun).toEqual([])
        const {node} = integrated.value;
        expect(node.root.findAllByType(Switch).length).toEqual(1);
        expect(node.root.findAllByProps({className:'reset-or-commit-inputs'}).length).toEqual(1)
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
        async function pushStepUI(){
          await act(async () => {
            bufferInputsContextData.value.setInputs(() => {
              return { 'bstep/b': [true], 'bstep/a': [true] }
            })
          });
        }
        describe('step marked pending', () => {
          const setStepComplete = awaitBefore(async () => {
            contextData.value.dispatch(setStatus(
              createStatusFixture(
                workflowHelpers.value.workflow,
                createStepStatusFixture({
                  name: 'bstep',
                  status: 'pending'
                })
              )
            ))
          })
          it('no validation errors, gets through TryCommitThenRun & triggers Run', async () => {
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
          it('with validation errors, no triggerRun, ', async () => {
            const {stateRef} = contextData.value;
            const {inputs} = bufferInputsContextData.value;
            
            await act(async () => {
              contextData.value.dispatch(
                addValidationErrors('bstep', 'the label', ['some error'])
              );
            });
            expect(stateRef.current.validationErrors).toEqual([["bstep","the label",["some error"]]]);

            // Triggered when setting inputs, with tryCommitThenRun also acted on within this await
            await pushStepUI()
            expect(stateRef.current.autoPassSteps).toEqual(['bstep']);
            expect(stateRef.current.tryCommitThenRun).toEqual([]);
            expect(stateRef.current.triggerRun).toEqual([]);
          })
        })
        describe('step completed', () => {
          const setStepComplete = awaitBefore(async () => {
            contextData.value.dispatch(setStatus(
              createStatusFixture(
                workflowHelpers.value.workflow,
                createStepStatusFixture({
                  name: 'bstep',
                  status: 'complete',
                  downloads: {output: 'https://download1'}
                })
              )
            ))
          })
          it('not triggered', async () => {
            const {stateRef} = contextData.value;
            const {inputs} = bufferInputsContextData.value;
            
            expect(stateRef.current.validationErrors).toEqual([]);
            
            // Triggered when setting inputs, with tryCommitThenRun also acted on within this await
            await pushStepUI()
            expect(stateRef.current.autoPassSteps).toEqual(['bstep']);
            expect(stateRef.current.tryCommitThenRun).toEqual([]);
            expect(stateRef.current.triggerRun).toEqual([]);
            const {node} = integrated.value;
            expect(node.root.findAllByProps({className:'reset-or-commit-inputs'}).length).toEqual(1)
          })
        })
      })
    })
    describe('for stepUserInputs w/out doc string', () => {
      
      stepName.replace(() => 'astep')
      
      it('not eligible - switch not shown', async () => {
        
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
});
