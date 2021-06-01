import React, {useContext} from 'react';
import {VulcanContext} from "./vulcan_context";
import {defaultWorkflow, Workflow} from "../api_types";

export function useWorkflow(): Workflow {
  const {state} = useContext(VulcanContext);
  return state.workflow || defaultWorkflow;
}
