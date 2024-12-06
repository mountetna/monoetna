import React, {useContext, useMemo} from 'react';
import {VulcanContext} from './vulcan_context';
import {defaultWorkflow, defaultWorkspace} from '../api_types';

export function useWorkspace() {
  const {state} = useContext(VulcanContext);
  const workflow = state.workflow || defaultWorkflow;
  const workspace = state.workspace || defaultWorkspace;
  const {stepsStatus} = state.status.steps;

  const hasPendingEdits = state.bufferedSteps.length > 0;

  const complete = useMemo(
    () => Object.keys(stepsStatus).length > 0 && Object.values(stepsStatus).every(({step}) => step.status === 'complete'),
    [stepsStatus]
  );

  return {workflow, workspace, hasPendingEdits, complete};
}
