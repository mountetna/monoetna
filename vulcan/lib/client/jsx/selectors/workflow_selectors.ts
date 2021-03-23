import {OUTPUT_COMPONENT, RUN, STATUS} from '../models/steps';
import {SessionStatusResponse, StepInput, StepStatus, Workflow, WorkflowStep} from "../api/types";
import {VulcanState} from "../reducers/vulcan_reducer";

export const workflowName = (workflow: Workflow | null | undefined) =>
    workflow && workflow.name ? workflow.name.replace('.cwl', '') : null;

export function workflowByName(name: string, state: VulcanState): Workflow | undefined {
    return state.workflows.find(w => workflowName(w) === name);
}

export function stepOfStatus(stepStatus: StepStatus, workflow: Workflow): WorkflowStep | undefined {
    return workflow.steps[0].find(s => s.name === stepStatus.name);
}

export function statusOfStep(step: WorkflowStep | string, status: VulcanState['status']): StepStatus | undefined {
    const stepName = (typeof step === 'string' ? step : step.name);
    return status[0].find(s => s.name === stepName);
}

export const uiOutputOfStep = (step: WorkflowStep) => {
    if(step.run && step.run.startsWith(RUN.UI_OUTPUT))
        return step.run.split('/')[1].replace('.cwl', '');
    return null;
}

export const uiQueryOfStep = (step: WorkflowStep) => {
    if(step.run && step.run.startsWith(RUN.UI_QUERY))
        return step.run.split('/')[1].replace('.cwl', '');
    return null;
}

export const uiStepInputDataLink = (step: WorkflowStep, status: VulcanState['status']): {[k: string]: string} => {
    // Pull out any previous step's output data link that is a required
    //   input into this UI step.
    const result: {[k: string]: string} = {};
    step.in.forEach(input => {
        const [stepName, outputKey] = input.source;
        const stepStatus = statusOfStep(stepName, status);
        if (stepStatus && stepStatus.downloads) {
            result[input.id] = stepStatus.downloads[outputKey];
        }
    })

    return result;
};

export const shouldDownload = (url: string, workflow: Workflow, step: WorkflowStep | undefined, status: SessionStatusResponse['status']) => {
    if (step == null) return false;
    return !!findDependentLinks(workflow, step, status).find(({ inputtingStep, outputtingStepStatus, outputName }) => {
        return !!(outputtingStepStatus.downloads && outputtingStepStatus.downloads[outputName] === url) &&
            isDataConsumer(inputtingStep);
    });
}

interface DependentLink {
    inputtingStep: WorkflowStep,
    stepInput: StepInput,
    outputtingStepStatus: StepStatus,
    outputName: string,
}

export const findDependentLinks = (workflow: Workflow, step: WorkflowStep, status: VulcanState['status']): DependentLink[] => {
    // Only download step data if its output is an input to
    //   a UI INPUT widget or a UI OUTPUT step that is not a LINK
    return workflow.steps[0].reduce((acc, otherStep) => {
        return acc.concat(dependentLinks(otherStep, step, status))
    }, [] as DependentLink[]);
};

export const findDependentSteps = (workflow: Workflow, step: WorkflowStep, status: VulcanState['status']) => {
    // Only download step data if its output is an input to
    //   a UI INPUT widget or a UI OUTPUT step that is not a LINK
    return workflow.steps[0].filter((s) => {
        return (dependentLinks(s, step, status).length > 0);
    });
};

export const dependentLinks = (step: WorkflowStep, dep: WorkflowStep, status: VulcanState['status']): DependentLink[] => {
    const result: DependentLink[] = [];
    const depStatus = statusOfStep(dep, status);

    if (depStatus) {
        step.in.forEach(stepInput => {
            const [inStepName, inOutputName] = stepInput.source;
            if (inStepName !== dep.name) return;

            dep.out.forEach(output => {
                if (inOutputName === output) {
                    result.push({
                        outputName: output,
                        outputtingStepStatus: depStatus,
                        inputtingStep: step,
                        stepInput,
                    });
                }
            })
        });
    }

    return result;
};

export const isDataConsumer = (step: WorkflowStep) => uiQueryOfStep(step) != null ||
    ([null, OUTPUT_COMPONENT.LINK].indexOf(uiOutputOfStep(step)) === -1)

export const uiStepInputDataRaw = (step: WorkflowStep, status: VulcanState['status'], data: VulcanState['data']): {[k: string]: any} => {
    // Pull out any previous step's output data link that is a required
    //   input into this UI step.
    const result: {[k: string]: any} = {};

    const urls = uiStepInputDataLink(step, status);
    Object.keys(urls).forEach(id => {
        const url = urls[id];
        const d = data[url];
        if (d) {
            result[id] = d;
        }
    })

    return result;
};

export const inputNameOfReference = (ref: [string] | [string, string]) => {
    if (ref.length === 1) return ref[0];
    return ref.join('/');
}

export const inputNamesOfStepOutputs = (step: WorkflowStep) => {
    return step.out.map(name => inputNameOfReference([step.name, name]))
}

export const defaultInputValues = (workflow: Workflow) => {
    return Object.keys(workflow.inputs).reduce((result, inputName) => {
        if (null != workflow.inputs[inputName].default) {
            result[inputName] = workflow.inputs[inputName].default;
        }
        return result;
    }, {});
};


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
        uiStepInputDataLink({step: workflowStep, status, pathIndex})
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

