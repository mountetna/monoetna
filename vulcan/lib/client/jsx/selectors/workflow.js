import {STATUS} from '../models/steps';
import {
  hasUiInput,
  hasUiOutput,
  workflowName as workflowNameShortener,
  uiInputDataReady
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

export const nextUiStepsSelector = (context) => {
  let {status, workflow, pathIndex} = context;
  return status[pathIndex]
    .map((s, index) => {
      let workflowStep = workflow.steps[pathIndex][index];
      if (
        STATUS.PENDING === s.status &&
        hasUiInput(workflowStep) &&
        uiInputDataReady({step: workflowStep, status, pathIndex})
      )
        return {...workflowStep, index};
    })
    .filter((s) => s);
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
