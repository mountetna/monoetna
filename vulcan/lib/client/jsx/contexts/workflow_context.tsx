import React, {useContext, useMemo} from 'react';
import {VulcanContext} from "./vulcan_context";
import {defaultWorkflow, Workflow} from "../api_types";
import {allWorkflowPrimaryInputSources, statusOfStep, uiOutputOfStep} from "../selectors/workflow_selectors";

export function useWorkflow() {
  const {state} = useContext(VulcanContext);
  const workflow = state.workflow || defaultWorkflow;
  const {status, session} = state;
  const {inputs} = session;
  const {steps} = workflow;

  const hasPendingEdits = state.bufferedSteps.length > 0;

  // We are done once every step either has a download or that step is a uiOutput.
  const complete = useMemo(
    () =>
      steps[0].every(
        (step) => uiOutputOfStep(step) || statusOfStep(step, status)?.downloads
      ),
    [steps, status]
  );

  return {workflow, hasPendingEdits, complete};
}
