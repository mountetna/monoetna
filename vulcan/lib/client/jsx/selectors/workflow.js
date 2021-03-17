import {STATUS} from '../models/steps';
import {
  hasUiInput,
  hasUiOutput,
  workflowName as workflowNameShortener
} from '../utils/workflow';

export const completedUiStepsSelector = (context) => {
  let {status, workflow, pathIndex} = context;
  return status[pathIndex]
    .map((step, index) => {
      let workflowStep = workflow.steps[pathIndex][index];
      if (STATUS.COMPLETE === step.status && hasUiInput(workflowStep)) {
        return {
          step: workflowStep,
          index
        };
      }
    })
    .filter((s) => s);
};

export const completedUiOutputsSelector = (context) => {
  let {status, workflow, pathIndex} = context;
  return status[pathIndex]
    .map((step, index) => {
      let workflowStep = workflow.steps[pathIndex][index];
      if (STATUS.COMPLETE === step.status && hasUiOutput(workflowStep)) {
        return {
          step: workflowStep,
          index
        };
      }
    })
    .filter((s) => s);
};

export const nextUiStepIndexSelector = (context) => {
  let {status, workflow, pathIndex} = context;
  return status[pathIndex].findIndex((s, index) => {
    let workflowStep = workflow.steps[pathIndex][index];
    return STATUS.PENDING === s.status && hasUiInput(workflowStep);
  });
};

export const errorStepsSelector = (context) => {
  let {status, workflow, pathIndex} = context;
  return status[pathIndex]
    .map((step, index) => {
      let workflowStep = workflow.steps[pathIndex][index];
      if (STATUS.ERROR === step.status) {
        return {
          step: workflowStep,
          index
        };
      }
    })
    .filter((s) => s);
};

export const workflowByName = ({workflows, workflowName}) =>
  workflows.find((w) => workflowNameShortener(w) === workflowName);
