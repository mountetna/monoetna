import React, {useContext, useMemo} from 'react';
import {VulcanContext} from "./vulcan_context";
import {defaultWorkflow} from "../api_types";

export function useWorkflow() {
  const {state} = useContext(VulcanContext);
  const workflow = state.workflow || defaultWorkflow;
  const {status} = state;

  const hasPendingEdits = state.bufferedSteps.length > 0;

  const complete = useMemo(
    () => status[0].length > 0 && status[0].every(({status}) => status === 'complete'),
    [status]
  );

  return {workflow, hasPendingEdits, complete};
}
