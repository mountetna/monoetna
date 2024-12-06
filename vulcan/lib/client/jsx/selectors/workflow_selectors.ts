import * as _ from 'lodash';
import {
  OUTPUT_COMPONENT,
  RUN,
  SessionStatusResponse,
  STATUS,
  StepInput,
  StepStatus,
  VulcanFigure,
  VulcanFigureSession,
  VulcanStorage,
  Workflow,
  WorkflowsResponse,
  WorkflowInput,
  WorkspaceStep,
  Workspace,
  VulcanConfigElement,
  defaultWorkflow
} from '../api_types';
import {VulcanState} from '../reducers/vulcan_reducer';
import {
  DataEnvelope,
  WorkflowStepGroup
} from '../components/workflow/user_interactions/inputs/input_types';
import {useMemo} from 'react';
import {
  mapSome,
  Maybe,
  maybeOfNullable,
  withDefault
} from 'etna-js/selectors/maybe';

export const workflowName = (workflow: Workflow | null | undefined) =>
  workflow && workflow.name ? workflow.name : null;

export function workflowByName(
  name: string,
  state: VulcanState
): Workflow | undefined {
  return state.workflows.find((w) => workflowName(w) === name);
}

export const workflowId = (workflow: Workflow | null | undefined) =>
  workflow && workflow.id ? workflow.id : null;

export function workflowById(
  id: any,
  state: VulcanState
): Workflow | undefined {
  return state.workflows.find((w) => workflowId(w) === id);
}

export function workflowByIdFromWorkflows(
  id: Workspace['workflow_id'],
  workflows: WorkflowsResponse
): Workflow | undefined {
  return id ? workflows.find((w) => workflowId(w) === id) : defaultWorkflow;
}

export const workspaceId = (workspace: Workspace | null | undefined) =>
  workspace && workspace.workspace_id ? workspace.workspace_id : null

export function workspaceById(
  id: any,
  state: VulcanState
): Workspace | undefined {
  return state.workspaces.find((w) => workspaceId(w) === id);
}

export function allUISteps(
  state: VulcanState
): string[] {
  return state.workspace ? Object.keys(state.workspace.vulcan_config) : [];
}

export function configUISteps(
  state: VulcanState
): string[] {
  if (!state.workspace) return []
  let uiSteps: string[] = []
  for (const [key,value] of Object.entries(state.workspace.vulcan_config)) {
    // Has params output
    // (Param outputs are to be exclusively adjusted by UIs for config elements.)
    if (value.output?.params) uiSteps.push(key)
  }
  return uiSteps
}

export function inputUISteps(
  state: VulcanState
): string[] {
  if (!state.workspace) return []
  let uiSteps: string[] = []
  for (const [key,value] of Object.entries(state.workspace.vulcan_config)) {
    // No params output, but has files output.  (MIGHT need to change this in the future)
    if (!value.output?.params && value.output?.files) uiSteps.push(key)
  }
  return uiSteps
}

// MIGHT need to change this in the future for e.g. optional usage of a plot output by drawing a box to select observations for a followup.
export function outputUISteps(
  state: VulcanState
): string[] {
  if (!state.workspace) return []
  let uiSteps: string[] = []
  for (const [key,value] of Object.entries(state.workspace.vulcan_config)) {
    // No outputs.  (MIGHT need to change this in the future)
    if (!value.output?.params && !value.output?.files) uiSteps.push(key)
  }
  return uiSteps
}

export function allUIFileInputs(
  workspace: Workspace
): string[] {
  if (!workspace) return [];
  let files: string[] = [];
  for (const value of Object.values(workspace.vulcan_config)) {
    if (value.input?.files) files.concat(value.input.files)
  };
  return files;
}

export function shouldDownloadOutput(
  workspace: Workspace,
  outputName: string
) {
  // Should download if output is an input to a UIstep -- which we can know from the vulcan config
  return allUIFileInputs(workspace).includes(outputName)
}

// export function stepOfStatus(
//   stepStatus: StepStatus | string,
//   workspace: Workspace
// ): WorkflowStep | undefined {
//   const stepName =
//     typeof stepStatus === 'string' ? stepStatus : stepStatus.name;
//   return workspace.steps[0].find((s) => s.name === stepName);
// }

export function statusOfStep(
  step: string,
  status: VulcanState['status']
): StepStatus | undefined {
  // const stepName = typeof step === 'string' ? step : step.name;
  const stepName = step;
  return status[0] ? status[0].find((s) => s.name === stepName) : undefined;
}

export function statusStringOfStepOrGroupedStep(
  step: WorkspaceStep | WorkflowStepGroup,
  workflow: Workflow,
  status: VulcanState['status']
) {
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

export function labelOfStepOrGroupedStep(
  step: WorkspaceStep | WorkflowStepGroup
) {
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

export const uiOutputOfStep = (step: WorkspaceStep) => {
  if (step.run && step.run.startsWith(RUN.UI_OUTPUT)) {
    return step.run.split('/')[1].replace('.cwl', '');
  }
  return null;
};

export const uiQueryOfStep = (step: WorkspaceStep) => {
  if (step.run && step.run.startsWith(RUN.UI_QUERY)) {
    return step.run.split('/')[1].replace('.cwl', '');
  }
  return null;
};

export const stepInputDataUrls = (
  step: WorkspaceStep,
  status: VulcanState['status']
): {[k: string]: string} => {
  // Pull out any previous step's output data link that is a required
  //   input into this UI step.
  const result: {[k: string]: string} = {};
  step.in.forEach((input) => {
    const [stepName, outputKey] = splitSource(input.source);

    if (stepName != null) {
      const stepStatus = statusOfStep(stepName, status);

      if (stepStatus && stepStatus.downloads) {
        result[input.id] = stepStatus.downloads[outputKey];
      }
    }
  });

  return result;
};

export function allWorkspacePrimaryInputSources(workspace: Workspace): string[] {
  let params = [] as string[];
  // ONLY ui for config elements have param outputs.
  Object.values(workspace.vulcan_config).map(val => {
    if (val.output && val.output.params) {
      params = params.concat(val.output.params);
    }
  })
  return params;
}

export const sourceNameOfReference = (ref: [string, string]) => {
  return ref.join('/');
};

export function inputSourcesOfStep(step: WorkspaceStep) {
  return step.out.map((outputName) =>
    sourceNameOfReference([step.name, outputName])
  );
}

export function inputSourcesOfPrimaryInputs(workflow: Workflow) {
  return Object.keys(workflow.inputs);
}

export function allWorkflowInputSources(workflow: Workflow): string[] {
  return allWorkflowPrimaryInputSources(workflow).concat(
    workflow.steps[0].reduce((acc, step) => {
      if (uiQueryOfStep(step)) {
        acc.push(...inputSourcesOfStep(step));
      }

      return acc;
    }, [] as string[])
  );
}

export function inputValueNonEmpty(
  val: Maybe<any>,
  disallow_empty_string = false,
  disallow_empty_array = true
): boolean {
  return withDefault(
    mapSome(
      val,
      (inner: any) =>
        (typeof inner !== 'number' || !isNaN(inner)) &&
        inner != null &&
        !_.isEqual(inner, ['']) &&
        (!disallow_empty_array || !_.isEqual(inner, [])) &&
        (!disallow_empty_string || !_.isEqual(inner, ''))
    ),
    false
  );
}

export function allDataNonEmpty(data: ([any] | null)[]) {
  return data.every((v) => v && inputValueNonEmpty(v[0]));
}

export function isNullish(value: any) {
  return value == null || value === 'null' || value === '';
}

export function dataOfSource(
  source: string,
  workflow: Workflow | null,
  status: VulcanState['status'],
  data: VulcanState['data'],
  session: VulcanState['session']
): [any] | null {
  if (!workflow) return null;
  const [originalStepName, outputName] = splitSource(source);
  const originalStep = originalStepName
    ? stepOfStatus(originalStepName, workflow)
    : null;
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

export function allExpectedOutputSources(
  step: VulcanConfigElement | WorkflowStepGroup
): string[] {
  if ('steps' in step) {
    return step.steps
      .map(allExpectedOutputSources)
      .reduce((a, b) => [...a, ...b], []);
  } else {
    if (!step.output) return [];
    let outputs: string[] = [];
    if (step.output.params) {
      outputs = outputs.concat(step.output.params);
    }
    if (step.output.files) {
      outputs = outputs.concat(step.output.files);
    }
    return outputs;
  }
}

export function allSourcesForStepName(
  name: string | null,
  workspace: Workspace | null
): string[] {
  if (!workspace) return [];
  if (!name) return allWorkspacePrimaryInputSources(workspace);
  if (!Object.keys(workspace.vulcan_config).includes(name)) return [];
  const step = workspace.vulcan_config[name];
  return allExpectedOutputSources(step).map(val => sourceNameOfReference([name, val]));
}

export const isDataConsumer = (step: WorkspaceStep) =>
  uiQueryOfStep(step) != null ||
  [null, OUTPUT_COMPONENT.LINK].indexOf(uiOutputOfStep(step)) === -1;

export function shouldDownloadStep(
  stepName: string,
  workflow: Workflow,
  outputName: string
) {
  return !!workflow.steps[0].find((step) => {
    if (!isDataConsumer(step)) return false;
    return !!step.in.find(({source}) => {
      const [sourceStep, sourceOutputName] = splitSource(source);
      return sourceStep === stepName && sourceOutputName == outputName;
    });
  });
}

export const stepInputDataRaw = (
  step: WorkspaceStep,
  status: VulcanState['status'],
  data: VulcanState['data'],
  session: VulcanState['session']
): {[k: string]: any} => {
  // Pull out any previous step's output data link that is a required
  //   input into this UI step.
  const result: {[k: string]: any} = {};

  const urls = stepInputDataUrls(step, status);
  step.in.forEach((input) => {
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
  });

  return result;
};

export function isPendingUiQuery(
  step: WorkspaceStep,
  status: VulcanState['status'],
  data: VulcanState['data'],
  session: VulcanState['session']
) {
  const bufferedData = stepInputDataRaw(step, status, data, session);
  return (
    uiQueryOfStep(step) &&
    statusOfStep(step, status)?.status == 'pending' &&
    step.in.every(({id}) => id in bufferedData)
  );
}

export const stepOutputs = (step: WorkspaceStep | undefined | '') => {
  if (!step) return [];

  return step.out.map((name) => sourceNameOfReference([step.name, name]));
};

export const sourceNamesOfStep = (step: WorkspaceStep) => {
  // return stepOutputs(step); Originally this...

  if (step.out.length === 0) return [];

  // Assign one output per CWL step so that any downloaded input data is grabbed.
  // But only take one output so that the input widget itself is not
  //   rendered multiple times in the UI.
  return [sourceNameOfReference([step.name, step.out[0]])];
};

export function missingOutputsForStep(
  step: WorkspaceStep,
  inputs: VulcanState['session']['inputs']
): {[k: string]: null} {
  const result: {[k: string]: null} = {};

  step.out.forEach((outputName) => {
    const source = sourceNameOfReference([step.name, outputName]);
    if (source in inputs) {
      return;
    }

    result[source] = null;
  });

  return result;
}

export function completedSteps(
  workflow: Workflow,
  status: VulcanState['status']
): WorkspaceStep[] {
  return workflow.steps[0].filter(
    (step) => statusOfStep(step, status)?.status === 'complete'
  );
}

export function completedUiOutputSteps(
  workflow: Workflow,
  status: VulcanState['status']
): WorkspaceStep[] {
  return completedSteps(workflow, status).filter(
    (step) => !!uiOutputOfStep(step)
  );
}

export function erroredSteps(
  workflow: Workflow,
  status: VulcanState['status']
) {
  return workflow.steps[0]
    .map((step, index) => ({step, index}))
    .filter(({step}) => statusOfStep(step, status)?.status === 'error');
}

export function pendingSteps(
  workflow: Workflow,
  status: VulcanState['status']
): WorkspaceStep[] {
  return workflow.steps[0].filter(
    (step) => statusOfStep(step, status)?.status === 'pending'
  );
}

export function hasNoRunningSteps(status: VulcanState['status']): boolean {
  return status[0].every((s) => s.status !== 'running');
}

export const inputGroupName = (name: string) => {
  let groupName = name.split('__')[0];
  if (groupName === name) return null;

  return groupName;
};

export function groupUiSteps(uiSteps: WorkspaceStep[]): WorkflowStepGroup[] {
  const map: {[k: string]: WorkflowStepGroup} = {};
  const result: WorkflowStepGroup[] = [];

  uiSteps.forEach((step) => {
    if ('isGroup' in step) {
      throw new Error('Cannot group a grouped input!' + JSON.stringify(step));
    }

    const groupName = inputGroupName(step.name);

    if (groupName == null) {
      result.push(
        (map[step.name] = {label: step.label || step.name, steps: [step]})
      );
      return;
    }

    if (!(groupName in map)) {
      result.push((map[groupName] = {label: groupName, steps: []}));
    }

    const group = map[groupName];
    group.steps.push(step);
  });

  return result;
}

export function filterEmptyValues(values: {[k: string]: any}): {
  [k: string]: any;
} {
  const result: {[k: string]: any} = {};

  Object.keys(values).forEach((k) => {
    const val = values[k];
    if (inputValueNonEmpty(maybeOfNullable(val))) result[k] = val;
  });

  return result;
}

export function useMemoized<P1, R>(f: (p1: P1) => R, p1: P1): R {
  return useMemo(() => f(p1), [f, p1]);
}

export function selectSession(
  figureResponse: VulcanFigureSession
): VulcanStorage {
  return {
    project_name: figureResponse.project_name,
    workflow_name: figureResponse.workflow_name,
    key: figureResponse.key,
    inputs: {...figureResponse.inputs},
    reference_figure_id: figureResponse.id
  };
}

export function selectFigure(
  figureResponse: VulcanFigureSession
): VulcanFigure {
  return {
    figure_id: figureResponse.figure_id,
    title: figureResponse.title,
    author: figureResponse.author,
    inputs: {...figureResponse.inputs},
    tags: [...(figureResponse.tags || [])],
    workflow_snapshot: figureResponse.workflow_snapshot,
    id: figureResponse.id
  };
}

export function mergeInputsWithDefaults(
  workflowInputs: Workflow['inputs'],
  sessionInputs: VulcanStorage['inputs'],
  currentInputs: VulcanStorage['inputs']
) {
  let withDefaults: DataEnvelope<Maybe<any>> = {};
  Object.keys(workflowInputs).forEach((inputName) => {
    if (!(inputName in sessionInputs) && !(inputName in currentInputs)) {
      withDefaults[inputName] = maybeOfNullable(
        workflowInputs[inputName].default
      );
    }
  });

  return withDefaults;
}

export function defaultInputs(workflow: Workflow) {
  return Object.entries(
    mergeInputsWithDefaults(workflow.inputs, {}, {})
  ).reduce(
    (
      acc: DataEnvelope<Maybe<any>>,
      [inputName, value]: [string, Maybe<any>]
    ) => {
      acc[inputName] = withDefault(value, null);

      return acc;
    },
    {}
  );
}
