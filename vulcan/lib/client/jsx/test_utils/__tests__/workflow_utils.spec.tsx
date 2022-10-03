import {awaitBefore, integrateElement, setupBefore} from '../integration';
import {useWorkflowUtils} from '../workflow_utils';
import {useContext} from 'react';
import {VulcanContext} from '../../contexts/vulcan_context';
import {stepInputDataRaw} from '../../selectors/workflow_selectors';

describe('workflow_utils', () => {
  describe('for a vulcan context', () => {
    const integrated = setupBefore(() => integrateElement());
    const workflowUtils = setupBefore(() => integrated.value.runHook(() => useWorkflowUtils()));
    const context = setupBefore(() => integrated.value.runHook(() => useContext(VulcanContext)));

    describe('for some steps flowing into another', () => {
      awaitBefore(async () => {
        workflowUtils.value.setWorkflow('test');
        workflowUtils.value.addStep('a', {
          out: ['out'],
        });
        workflowUtils.value.addStep('b', {
          out: ['out'],
        });
        workflowUtils.value.addStep('c', {
          in: [{ source: 'a/out', id: 'a' }, { source: 'b/out', id: 'b' }],
        });
      });

      describe('setting downloads for the inputs', () => {
        awaitBefore(async () => {
          workflowUtils.value.forceDownloadedData('a/out', 'avalue');
          workflowUtils.value.forceDownloadedData('b/out', 'bvalue');
        });

        it('prepares the downloads successfully as inputs', () => {
          const {status, data, session} = context.value.stateRef.current;
          expect(stepInputDataRaw(workflowUtils.value.steps['c'], status, data, session)).toEqual({
            b: 'bvalue',
            a: 'avalue',
          });
        });
      });
    });
    describe('workflow with no vignette', () => {
      awaitBefore(async () => {
        workflowUtils.value.setWorkflow('test');
        workflowUtils.value.addStep('a', {
          out: ['out'],
        });
      });
      const context = setupBefore(() => integrated.value.runHook(() => useContext(VulcanContext)));
      it('has no vignette', () => {
        expect(context.value.state.workflow?.vignette).toBeUndefined();
      });
      describe('setVignette', () => {
        const addVignette = awaitBefore(async () => {
          workflowUtils.value.setVignette('I am a helpful vignette.');
        });
        const context = setupBefore(() => integrated.value.runHook(() => useContext(VulcanContext)));
        it('adds a vignette to the current workflow', async () => {
          expect(context.value.state.workflow?.vignette).toEqual('I am a helpful vignette.');
        });
      });
    });
  });
});