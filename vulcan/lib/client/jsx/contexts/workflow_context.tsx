import React, {useContext, useMemo} from 'react';
import {VulcanContext} from "./vulcan_context";
import {defaultWorkflow, Workflow} from "../api_types";
import {allWorkflowPrimaryInputSources, statusOfStep, uiOutputOfStep} from "../selectors/workflow_selectors";

export function useWorkflow() {
  const {state} = useContext(VulcanContext);
  const workflow = state.workflow || defaultWorkflow;
  const {status, session} = state;

  const hasPendingEdits = state.bufferedSteps.length > 0;

  const complete = useMemo(
    () => status[0].every(({status}) => status === 'complete'),
    [status]
  );

  return {workflow, hasPendingEdits, complete};
}
