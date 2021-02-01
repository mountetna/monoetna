export const FETCH_WORKFLOWS = 'FETCH_WORKFLOWS';
export const FETCH_WORKFLOW = 'FETCH_WORKFLOW';
export const SUBMIT_INPUTS = 'SUBMIT_INPUTS';

export const fetchWorkflows = (workflows) => {
  return {
    type: FETCH_WORKFLOWS,
    workflows: workflows
  };
};

export const fetchWorkflow = (workflow) => {
  return {
    type: FETCH_WORKFLOW,
    workflow: workflow
  };
};

export const submitInput = (status) => {
  return {
    type: SUBMIT_INPUTS,
    workflow: status
  };
};
