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
    return !!findDependentLinks(workflow, step, status).find(({ inputtingStep, source }) => {
        const { outputtingStepStatus, outputName } = source;
        return !!(outputtingStepStatus.downloads && outputtingStepStatus.downloads[outputName] === url) &&
            isDataConsumer(inputtingStep);
    });
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

type IOLink = StepToStepLink | PrimaryInputToStepLink;

export const findDependentLinks = (workflow: Workflow, step: WorkflowStep, status: VulcanState['status']): StepToStepLink[] => {
    return workflow.steps[0].reduce((acc, otherStep) => {
        return acc.concat(linksBetweenSteps(otherStep, step, status))
    }, [] as StepToStepLink[]);
};

export function findDependencyLinks(workflow: Workflow, step: WorkflowStep, status: VulcanState['status']): IOLink[] {
    return step.in.reduce((acc, stepIn) => {
        const link = ioLinkOf(step, stepIn.id, status, workflow)
        if (link) acc.push(link);
        return acc;
    }, [] as IOLink[]);
}

export function findPrimaryInputDependentLinks(workflow: Workflow): PrimaryInputToStepLink[] {
    return workflow.steps[0].reduce((links, step) => {
        step.in.forEach(stepIn => {
            const link = ioLinkOf(step, stepIn.id, defaultVulcanState.status, workflow);
            if (link) {
                links.push(link);
            }
        });

        return links;
    }, [] as PrimaryInputToStepLink[]);
}

export function findFulfilledSteps(workflow: Workflow, step: WorkflowStep, status: VulcanState['status']) {
    const result: {[k: string]: [StepStatus, IOLink[]]} = {};

    workflow.steps[0].filter(step => {
        step.in.forEach(stepInput => {
            // const {outputtingStepStatus} = dependentLink()
        })
    })

    return result;
}

export function ioLinkOf(inStep: WorkflowStep, inputName: string,
                         status: VulcanState['status'], depSrc: Workflow | WorkflowStep): IOLink | null {
    const stepInput = inStep.in.find(({id}) => id === inputName);
    if (stepInput == null) return null;
    const [depName, depOutput] = stepInput.source;

    if (depName !== "primary_inputs") {
        const depStatus = statusOfStep(depName, status);
        const depStep = 'steps' in depSrc ? stepOfStatus(depName, depSrc) : depSrc;

        if (!depStatus || !depStep) return null;
        if (depStep.out.indexOf(depOutput) === -1) return null;

        return {
            source: {
                outputName: depOutput,
                outputtingStepStatus: depStatus,
            },
            inputtingStep: inStep,
            stepInput,
        }
    } else {
        return {
            // primary input source names do not include a slash
            source: depOutput,
            inputtingStep: inStep,
            stepInput,
        }
    }
}

export const linksBetweenSteps = (inStep: WorkflowStep, outStep: WorkflowStep, status: VulcanState['status']): StepToStepLink[] => {
    const result: StepToStepLink[] = [];
    const depStatus = statusOfStep(outStep, status);

    if (depStatus) {
        inStep.in.forEach(stepInput => {
            const link = ioLinkOf(inStep, stepInput.id, status, outStep);
            if (link != null) result.push(link as StepToStepLink);
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

