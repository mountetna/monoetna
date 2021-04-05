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
  if (step.run && step.run.startsWith(RUN.UI_OUTPUT))
    return step.run.split('/')[1].replace('.cwl', '');
  return null;
}

export const uiQueryOfStep = (step: WorkflowStep) => {
  if (step.run && step.run.startsWith(RUN.UI_QUERY))
    return step.run.split('/')[1].replace('.cwl', '');
  return null;
}

export const stepInputDataUrls = (step: WorkflowStep, status: VulcanState['status']): { [k: string]: string } => {
  // Pull out any previous step's output data link that is a required
  //   input into this UI step.
  const result: { [k: string]: string } = {};
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

export function dependentStepConsumersOf(startingSource: string, workflow: Workflow, shallow = false) {
  const result: WorkflowStep[] = [];
  const dependents = workflow.dependencies_of_outputs[startingSource];

  if (dependents != null) {
    result.push(...workflow.steps[0].filter(s => !!uiOutputOfStep(s)).filter(
        s => !!s.in.find(({source}) => source === startingSource || dependents.indexOf(source) !== -1)));
  }

  dependents.map(s => stepOfSource(s)).forEach(stepName => {
    if (!stepName) return;
    const step = stepOfStatus(stepName, workflow);

    if (step) {
      if (result.indexOf(step) === -1) result.push(step);
    }
  });

  if (shallow) {
    return result.filter(step => {
      return step.in.find(({source}) => stepOfSource(source) == stepOfSource(startingSource));
    })
  }

  return result;
}

export const shouldDownload = (url: string, workflow: Workflow, step: WorkflowStep | undefined, status: SessionStatusResponse['status']) => {
  if (step == null) return false;
  const stepStatus = statusOfStep(step, status);
  if (stepStatus == null) return false;
  const {downloads} = stepStatus;
  if (downloads == null) return false;

  return !!step.out.find(outName => {
    const source = sourceNameOfReference([step.name, outName]);
    if (downloads[outName] !== url) return false;
    const deps = dependentStepConsumersOf(source, workflow, true);
    return !!deps.find(dep => isDataConsumer(dep));
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

export const stepInputDataRaw = (
    step: WorkflowStep,
    status: VulcanState['status'], data: VulcanState['data'], session: VulcanState['session'],
): { [k: string]: any } => {
  // Pull out any previous step's output data link that is a required
  //   input into this UI step.
  const result: { [k: string]: any } = {};

  const urls = stepInputDataUrls(step, status);
  step.in.forEach(input => {
    if (input.source in session.inputs) {
      result[input.id] = session.inputs[input.source];
    } else {
      if (input.id in urls) {
        const url = urls[input.id];
        if (url in data) {
          result[input.id] = data[urls[input.id]];
        }
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
  }, {} as { [k: string]: any });
};

export function missingUiQueryOutputs(step: WorkflowStep, inputs: VulcanState['inputs']): { [k: string]: null } {
  const result: { [k: string]: null } = {};

  step.out.forEach(outputName => {
    const source = sourceNameOfReference([step.name, outputName]);
    if (source in inputs) {
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
  return workflow.steps[0].map((step, index) => ({step, index}))
      .filter(({step}) => statusOfStep(step, status)?.status === 'complete');
}

export function erroredSteps(workflow: Workflow, status: VulcanState['status']) {
  return workflow.steps[0].map((step, index) => ({step, index}))
      .filter(({step}) => statusOfStep(step, status)?.status === 'error');
}


export function pendingSteps(workflow: Workflow, status: VulcanState['status']): UIStep[] {
  return workflow.steps[0].map((step, index) => ({step, index}))
      .filter(({step}) => statusOfStep(step, status)?.status === 'pending');
}

export const inputGroupName = (name: string) => {
  let groupName = name.split('__')[0];
  if (groupName === name) groupName = 'Inputs';

  groupName = groupName.replace(/_/g, ' ');

  return groupName;
};

export function groupUiSteps(uiSteps: UIStep[]): UIStep[] {
  const result: { [k: string]: UIStep } = {};

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

export function filterEmptyValues(values: { [k: string]: any }): { [k: string]: any } {
  const result: { [k: string]: any } = {};

  Object.keys(values).forEach(k => {
    const val = values[k];
    if (inputValueNonEmpty(val)) result[k] = val;
  })

  return result;
}

// For each key in changes, remove any dependencies_of_outputs for that key (dependents)
// in the existing object (copy on change), iff that dependent is not also in changes.
export function unsetDependentInputs(
    changes: { [k: string]: any },
    existing: { [k: string]: any },
    workflow: Workflow | null
): { [k: string]: any } {
  let result = existing;
  if (!workflow) return result;

  Object.keys(changes).forEach(changedSource => {
    const dependents = workflow.dependencies_of_outputs[changedSource];
    if (!dependents) return;

    dependents.forEach(dependentSource => {
      // In the case of forward change propagation, allow specifying two inputs at once that happen to have a dependency
      if (dependentSource in changes) return;
      result = {...result};
      delete result[dependentSource];
    });
  });

  return result;
}

export function findSourceDependencies(source: string, workflow: Workflow | null): string[] {
  if (!workflow) return [];
  return Object.keys(workflow.dependencies_of_outputs).filter(k => workflow.dependencies_of_outputs[k].includes(source));
}
