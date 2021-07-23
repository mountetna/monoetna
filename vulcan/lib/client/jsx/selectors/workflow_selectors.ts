import * as _ from 'lodash';
import {
  OUTPUT_COMPONENT, RUN, SessionStatusResponse, STATUS, StepInput, StepStatus, Workflow, WorkflowInput, WorkflowStep,
} from "../api_types";
import {VulcanState} from "../reducers/vulcan_reducer";
import {
   WorkflowStepGroup
} from "../components/workflow/user_interactions/inputs/input_types";
import {useMemo} from "react";
import {mapSome, Maybe, maybeOfNullable, withDefault} from "./maybe";

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
  return status[0] ? status[0].find(s => s.name === stepName) : undefined;
}

export function statusStringOfStepOrGroupedStep(step: WorkflowStep | WorkflowStepGroup, workflow: Workflow, status: VulcanState['status']) {
  if ('steps' in step) {
    let statusStr = null as string | null;
    for (let innerStep of step.steps) {
      let stepStatus = statusOfStep(innerStep.name, status);

      if (!stepStatus) {
        return STATUS.PENDING;
      }

      if (statusStr == null) statusStr = stepStatus.status;
      else if (statusStr !== stepStatus.status) return STATUS.PENDING;
    }

    return statusStr || STATUS.COMPLETE;
  }

  return statusOfStep(step, status)?.status || STATUS.PENDING;
}

export function labelOfStepOrGroupedStep(step: WorkflowStep | WorkflowStepGroup) {
  if ('steps' in step) {
    return step.label;
  }

  return step.label || step.name;
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

export function inputSourcesOfStep(step: WorkflowStep) {
  return step.out.map(outputName => sourceNameOfReference([step.name, outputName]));
}

export function inputSourcesOfPrimaryInputs(workflow: Workflow) {
  return Object.keys(workflow.inputs);
}

export function allWorkflowInputSources(workflow: Workflow): string[] {
  return allWorkflowPrimaryInputSources(workflow).concat(workflow.steps[0].reduce((acc, step) => {
    if (uiQueryOfStep(step)) {
      acc.push(...inputSourcesOfStep(step));
    }

    return acc;
  }, [] as string[]));
}

export function inputValueNonEmpty(val: Maybe<any>, disallow_empty_string = false): boolean {
  return withDefault(mapSome(val, inner =>
    (typeof inner !== "number" || !isNaN(inner)) &&
      inner != null &&
      !_.isEqual(inner, ['']) &&
      !_.isEqual(inner, []) &&
      (!disallow_empty_string || !_.isEqual(inner, ''))
    ), false);
}

export function allDataNonEmpty(data: ([any] | null)[]) {
  return data.every(v => v && inputValueNonEmpty(v[0]));
}

export function dataOfSource(source: string, workflow: Workflow | null, status: VulcanState['status'], data: VulcanState['data'], session: VulcanState['session']): [any] | null {
  if (!workflow) return null;
  const [originalStepName, outputName] = splitSource(source);
  const originalStep = originalStepName ? stepOfStatus(originalStepName, workflow) : null;
  if (!originalStep) return null;

  if (source in session.inputs) return [session.inputs[source]];
  const stepStatus = statusOfStep(originalStep, status);
  if (!stepStatus) return null;
  const {downloads} = stepStatus;
  if (!downloads) return null;
  const url = downloads[outputName];
  if (url in data) return data[url];
  return null;
}

export function allExpectedOutputSources(step: WorkflowStep | WorkflowStepGroup): string[] {
  if ('steps' in step) {
    return step.steps.map(allExpectedOutputSources).reduce((a, b) => [...a, ...b], []);
  } else {
    return step.out.map(outputName => sourceNameOfReference([step.name, outputName]));
  }
}

export function allSourcesForStepName(name: string | null, workflow: Workflow | null): string[] {
  if (!workflow) return [];
  if (!name) return allWorkflowPrimaryInputSources(workflow);
  const step = stepOfStatus(name, workflow);
  if (!step) return [];
  return allExpectedOutputSources(step);
}

export function shouldDownloadStep(stepName: string, workflow: Workflow, outputName: string) {
  return !!workflow.steps[0].find(step => {
    if (!isDataConsumer(step)) return false;
    return !!step.in.find(({source}) => {
      const [sourceStep, sourceOutputName] = splitSource(source);
      return sourceStep === stepName && sourceOutputName == outputName;
    });
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
    status: VulcanState['status'],
    data: VulcanState['data'],
    session: VulcanState['session'],
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

export function missingOutputsForStep(step: WorkflowStep, inputs: VulcanState['session']['inputs']): { [k: string]: null } {
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

export function completedUiOutputSteps(workflow: Workflow, status: VulcanState['status']): WorkflowStep[] {
  return completedSteps(workflow, status).filter((step) => !!uiOutputOfStep(step));
}

export function completedSteps(workflow: Workflow, status: VulcanState['status']): WorkflowStep[] {
  return workflow.steps[0].filter((step) => statusOfStep(step, status)?.status === 'complete');
}

export function erroredSteps(workflow: Workflow, status: VulcanState['status']) {
  return workflow.steps[0].map((step, index) => ({step, index}))
      .filter(({step}) => statusOfStep(step, status)?.status === 'error');
}


export function pendingSteps(workflow: Workflow, status: VulcanState['status']): WorkflowStep[] {
  return workflow.steps[0].filter(step => statusOfStep(step, status)?.status === 'pending');
}

export function hasNoRunningSteps(status: VulcanState['status']): boolean {
  return status[0].every(s => s.status !== "running");
}

export const inputGroupName = (name: string) => {
  let groupName = name.split('__')[0];
  if (groupName === name) return null;

  groupName = groupName.replace(/_/g, ' ');

  return groupName;
};

export function groupUiSteps(uiSteps: WorkflowStep[]): WorkflowStepGroup[] {
  const map: { [k: string]: WorkflowStepGroup } = {};
  const result: WorkflowStepGroup[] = [];

  uiSteps.forEach(step => {
    if ('isGroup' in step) {
      throw new Error('Cannot group a grouped input!' + JSON.stringify(step));
    }

    const groupName = inputGroupName(step.name);

    if (groupName == null) {
      result.push(map[step.name] = { label: "Inputs", steps: [step] });
      return;
    }

    if (!(groupName in map)) {
      result.push(map[groupName] = { label: groupName, steps: [] })
    }

    const group = map[groupName];
    group.steps.push(step);
  });

  return result;
}

export function filterEmptyValues(values: { [k: string]: any }): { [k: string]: any } {
  const result: { [k: string]: any } = {};

  Object.keys(values).forEach(k => {
    const val = values[k];
    if (inputValueNonEmpty(maybeOfNullable(val))) result[k] = val;
  })

  return result;
}

export function useMemoized<P1, R>(f: (p1: P1) => R, p1: P1): R {
  return useMemo(() => f(p1), [f, p1]);
}