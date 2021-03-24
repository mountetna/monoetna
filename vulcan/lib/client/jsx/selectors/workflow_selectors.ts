import * as _ from 'lodash';
import {OUTPUT_COMPONENT, RUN, STATUS, SessionStatusResponse, StepInput, StepStatus, Workflow, WorkflowStep} from "../api_types";
import {defaultVulcanState, VulcanState} from "../reducers/vulcan_reducer";

export const workflowName = (workflow: Workflow | null | undefined) =>
    workflow && workflow.name ? workflow.name.replace('.cwl', '') : null;

export function workflowByName(name: string, state: VulcanState): Workflow | undefined {
    return state.workflows.find(w => workflowName(w) === name);
}

export function stepOfStatus(stepStatus: StepStatus | string, workflow: Workflow): WorkflowStep | undefined {
    const stepName = (typeof stepStatus === 'string' ? stepStatus : stepStatus.name);
    return workflow.steps[0].find(s => s.name === stepName);
}

export function statusOfStep(step: WorkflowStep | string, status: VulcanState['status']): StepStatus | undefined {
    const stepName = (typeof step === 'string' ? step : step.name);
    return status[0].find(s => s.name === stepName);
}

const SOURCE_STR_REGEX = /^([^\/]+)\/[^\/]+$/;
export function stepOfSource(source: string): string | undefined {
    const match = source.match(SOURCE_STR_REGEX);
    if (match) return match[1];
    return undefined;
}

export function splitSource(source: string): [string | undefined, string] {
    const step = stepOfSource(source);
    if (step) source = source.slice(step.length + 1);
    return [step, source];
}


export const uiOutputOfStep = (step: WorkflowStep) => {
    if(step.run && step.run.startsWith(RUN.UI_OUTPUT))
        return step.run.split('/')[1].replace('.cwl', '');
    return null;
}

export const uiQueryOfStep = (step: WorkflowStep) => {
    if (step.run && step.run.startsWith(RUN.UI_QUERY))
        return step.run.split('/')[1].replace('.cwl', '');
    return null;
}

export const stepInputDataUrls = (step: WorkflowStep, status: VulcanState['status']): {[k: string]: string} => {
    // Pull out any previous step's output data link that is a required
    //   input into this UI step.
    const result: {[k: string]: string} = {};
    step.in.forEach(input => {
        const [stepName, outputKey] = splitSource(input.source);

        if (stepName != null) {
            const stepStatus = statusOfStep(stepName, status);

            if (stepStatus && stepStatus.downloads) {
                result[input.id] = stepStatus.downloads[outputKey];
            }
        }
    })

    return result;
};

export function allWorkflowPrimaryInputSources(workflow: Workflow): string[] {
    return Object.keys(workflow.inputs);
}

export function allWorkflowInputSources(workflow: Workflow): string[] {
    return allWorkflowPrimaryInputSources(workflow).concat(workflow.steps[0].reduce((acc, step) => {
        if (uiQueryOfStep(step)) {
            acc.push(...step.out.map(outputName => sourceNameOfReference([step.name, outputName])));
        }

        return acc;
    }, [] as string[]));
}

export function inputValueNonEmpty(val: any) {
    return val != null && !isNaN(val) && !_.isEqual(val, ['']) && !_.isEqual(val, []);
}

export const shouldDownload = (url: string, workflow: Workflow, step: WorkflowStep | undefined, status: SessionStatusResponse['status']) => {
    if (step == null) return false;
    const stepStatus = statusOfStep(step, status);
    if (stepStatus == null) return false;
    const { downloads } = stepStatus;
    if (downloads == null) return false;

    return !!sourceNamesOfStep(step).find(source => {
        if(downloads[source] !== url) return false;
        const deps = workflow.dependencies_of_outputs[source];
        if (!deps || deps.length === 0) return false;
        return !!deps.find(depSource => {
            const depStepName = stepOfSource(depSource);
            if (depStepName == null) return false;
            const depStep = stepOfStatus(depStepName, workflow);
            if (depStep == null) return false;
            return isDataConsumer(depStep);
        })
    })
}

interface StepToStepLink {
    inputtingStep: WorkflowStep,
    stepInput: StepInput,
    source: {
        outputtingStepStatus: StepStatus,
        outputName: string,
    };
}

interface PrimaryInputToStepLink {
    inputtingStep: WorkflowStep,
    stepInput: StepInput,
    source: string;
}

export const isDataConsumer = (step: WorkflowStep) => uiQueryOfStep(step) != null ||
    ([null, OUTPUT_COMPONENT.LINK].indexOf(uiOutputOfStep(step)) === -1)

export const stepInputDataRaw = (step: WorkflowStep, status: VulcanState['status'], data: VulcanState['data'], session: VulcanState['session']): {[k: string]: any} => {
    // Pull out any previous step's output data link that is a required
    //   input into this UI step.
    const result: {[k: string]: any} = {};

    const urls = stepInputDataUrls(step, status);
    step.in.forEach(input => {
        if (input.source in session.inputs) {
            result[input.id] = session.inputs[input.source];
        } else {
            if (input.id in urls) {
                result[input.id] = urls[input.id];
            }
        }
    })

    return result;
};

export const sourceNameOfReference = (ref: [string, string]) => {
    return ref.join('/');
}

export const sourceNamesOfStep = (step: WorkflowStep) => {
    return step.out.map(name => sourceNameOfReference([step.name, name]))
}

export const defaultInputValues = (workflow: Workflow) => {
    return Object.keys(workflow.inputs).reduce((result, inputName) => {
        if (null != workflow.inputs[inputName].default) {
            result[inputName] = workflow.inputs[inputName].default;
        }
        return result;
    }, {} as {[k: string]: any});
};


// export const completedUiStepsSelector = (context) => {
//   let {status, workflow, pathIndex} = context;
//   return status[pathIndex]
//     .map((step, index) => {
//       let workflowStep = workflow.steps[pathIndex][index];
//       if (STATUS.COMPLETE === step.status && hasUiInput(workflowStep)) {
//         return {
//           step: workflowStep,
//           index
//         };
//       }
//     })
//     .filter((s) => s);
// };
//
// export const completedUiOutputsSelector = (context) => {
//   let {status, workflow, pathIndex} = context;
//   return status[pathIndex]
//     .map((step, index) => {
//       let workflowStep = workflow.steps[pathIndex][index];
//       if (STATUS.COMPLETE === step.status && hasUiOutput(workflowStep)) {
//         return {
//           step: workflowStep,
//           index
//         };
//       }
//     })
//     .filter((s) => s);
// };
//
// export const nextUiStepsSelector = (context) => {
//   let {status, workflow, pathIndex} = context;
//   return status[pathIndex]
//     .map((s, index) => {
//       let workflowStep = workflow.steps[pathIndex][index];
//       if (
//         STATUS.PENDING === s.status &&
//         hasUiInput(workflowStep) &&
//         uiStepInputDataLink({step: workflowStep, status, pathIndex})
//       )
//         return {...workflowStep, index};
//     })
//     .filter((s) => s);
// };
//
// export const errorStepsSelector = (context) => {
//   let {status, workflow, pathIndex} = context;
//   return status[pathIndex]
//     .map((step, index) => {
//       let workflowStep = workflow.steps[pathIndex][index];
//       if (STATUS.ERROR === step.status) {
//         return {
//           step: workflowStep,
//           index
//         };
//       }
//     })
//     .filter((s) => s);
// };

