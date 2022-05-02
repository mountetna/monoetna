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
  setBufferedInput,
  setInputs,
  setSession
} from '../../actions/vulcan_actions';
import {
  BufferedInputsContext,
  WithBufferedInputs
} from '../input_state_management';
import {createSessionFixture} from '../../test_utils/fixtures';

describe('useDataBuffering', () => {
  beforeEach(() => {
    jest.useFakeTimers();
  });

  afterEach(() => {
    jest.useRealTimers();
  });

  const stepName = setupBefore(() => null as string | null);
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
  
  const setVignette = awaitBefore(async () => {
    workflowHelpers.value.setVignette('blah blah Primary inputs are skippable blah blah');
  }); 

  awaitBefore(async () => {
    expect(
      [...contextData.value.stateRef.current.bufferedSteps].sort()
    ).toEqual([stepName.value]);
  });

  describe('autoPass', () => {
    const doCancel = awaitBefore(async () => {
        bufferInputsContextData.value.cancelInputs();
      });  
    // describe('for primary inputs w/ vignette string', () => {
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
    describe('for stepUserInputs w/ vignette string', () => {
      // const addAutoPassableInput = awaitBefore(async () => {
      //   contextData.value.dispatch(setInputs({v: 2}));
      // });
      it('eligible - switch shown', () => {})
        expect(contextData.value.state.autoPassSteps).toEqual([])
        expect(contextData.value.state.tryCommitThenRun).toEqual([])
        expect(contextData.value.state.triggerRun).toEqual([])
        // expect(node.root.findAllByProps({className: 'auto-pass-toggle'}).length).toEqual(1);
      it('triggers TryCommitThenRun onBuffer if pending but not complete', () => {})
      it('clears TryCommitThenRun & triggers Run if commit is successful', () => {})
      it('clears TryCommitThenRun, if commit has validation errors', () => {})
    })
    describe('for stepUserInputs w/out vignette string', () => {
      stepName.replace(() => 'bstep')
      it('not eligible - switch not shown', () => {})
    })
  })
});
