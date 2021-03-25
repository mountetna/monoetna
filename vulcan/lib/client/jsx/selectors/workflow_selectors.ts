import * as _ from 'lodash';
import {
    OUTPUT_COMPONENT,
    RUN,
    SessionStatusResponse, StepInput,
    StepStatus,
    Workflow,
    WorkflowStep,
} from "../api_types";
import {VulcanState} from "../reducers/vulcan_reducer";
import {GroupedInputStep, UIStep} from "../components/workflow/user_interactions/inputs/input_types";

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

export function allDataNonEmpty(data: ([any] | null)[]) {
    return data.every(v => v && inputValueNonEmpty(v[0]));
}

export function dataOfSource(source: string, workflow: Workflow | null, status: VulcanState['status'], data: VulcanState['data'], session: VulcanState['session']): [any] | null {
    if (!workflow) return null;
    const [originalStepName, outputName] = splitSource(source);
    const originalStep = originalStepName ? stepOfStatus(originalStepName, workflow) : null;
    return originalStep ? [stepInputDataRaw(originalStep, status, data, session)[outputName]] : null;
}

export function allExpectedOutputSources(step: WorkflowStep | GroupedInputStep): string[] {
    if ('isGroup' in step) {
        return step.in.map(({source}) => source);
    } else {
        return step.out.map(outputName => sourceNameOfReference([step.name, outputName]));
    }
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

export const isDataConsumer = (step: WorkflowStep) => uiQueryOfStep(step) != null ||
    ([null, OUTPUT_COMPONENT.LINK].indexOf(uiOutputOfStep(step)) === -1)

export function isPendingUiQuery(step: WorkflowStep, status: VulcanState['status'], data: VulcanState['data'], session: VulcanState['session']) {
    const bufferedData = stepInputDataRaw(step, status, data, session);
    return uiQueryOfStep(step)
        && statusOfStep(step, status)?.status == 'pending'
        && step.in.every(({id}) => id in bufferedData)
}

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

export function missingUiQueryOutputs(step: WorkflowStep, inputs: VulcanState['inputs']): {[k: string]: null} {
    const result: {[k: string]: null} = {};

    step.out.forEach(outputName => {
        const source = sourceNameOfReference([step.name, outputName]);
        if (source in inputs)  {
            return;
        }

        result[source] = null;
    });

    return result;
}

export function completedUiOutputSteps(workflow: Workflow, status: VulcanState['status']): UIStep[] {
    return completedSteps(workflow, status).filter(({step}) => !!uiOutputOfStep(step));
}

export function completedSteps(workflow: Workflow, status: VulcanState['status']): UIStep[] {
    return workflow.steps[0].map((step, index) => ({ step, index }))
        .filter(({step}) => statusOfStep(step, status)?.status === 'complete');
}
export function erroredSteps(workflow: Workflow, status: VulcanState['status']) {
    return workflow.steps[0].map((step, index) => ({ step, index }))
        .filter(({step}) => statusOfStep(step, status)?.status === 'error');
}


export function pendingSteps(workflow: Workflow, status: VulcanState['status']): UIStep[] {
    return workflow.steps[0].map((step, index) => ({ step, index }))
        .filter(({step}) => statusOfStep(step, status)?.status === 'pending');
}

export const inputGroupName = (name: string) => {
  let groupName = name.split('__')[0];
  if (groupName === name) groupName = 'Inputs';

  groupName = groupName.replace(/_/g, ' ');

  return groupName;
};

export function groupUiSteps(uiSteps: UIStep[]): UIStep[] {
    const result: {[k: string]: UIStep} = {};

    uiSteps.forEach(uiStep => {
        if ('isGroup' in uiStep) {
            throw new Error('Cannot group a grouped input!' + JSON.stringify(uiStep));
        }

        const groupName = inputGroupName(uiStep.step.name);

        if (groupName === 'Inputs') {
            result[uiStep.step.name] = uiStep;
            return;
        }

        const group = result[groupName] = result[groupName] || {
            step: {
                name: groupName,
                isGroup: true,
                label: groupName,
                run: "", // group steps have no real run
                in: [],
            },
            index: uiStep.index,
        }

        group.step.in = group.step.in.concat(uiStep.step.out.map(name => ({
            id: sourceNameOfReference([uiStep.step.name, name]),
            source: sourceNameOfReference([uiStep.step.name, name]),
            doc: uiStep.step.doc,
            label: uiStep.step.label || sourceNameOfReference([uiStep.step.name, name]),
        })));
    });

    return Object.values(result);
}
