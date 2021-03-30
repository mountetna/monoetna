import {setStatus, setWorkflows} from "../actions/vulcan";
import {workflowsResponse} from "./fixtures/workflows-response";

export const setWorkflowsFixture = setWorkflows(workflowsResponse.workflows);
export const setStatusFixture = setStatus()
