import * as _ from 'lodash';
import {
  STATUS,
  StepStatus,
  Workflow,
  WorkflowsResponse,
  WorkspaceStep,
  Workspace,
  VulcanConfigElement,
  defaultWorkflow,
  WorkspaceStatus,
  StatusString,
  VulcanConfig,
  VulcanConfigRaw,
  FlatParams,
  MultiFileContentResponse,
  MultiFileContent,
  RunStatus,
  StatusStringBroaden,
  WorkspaceRaw
} from '../api_types';
import {VulcanState} from '../reducers/vulcan_reducer';
import {
  WorkspaceStepGroup
} from '../components/workflow/user_interactions/inputs/input_types';
import {useMemo} from 'react';
import {
  isSome,
  mapSome,
  Maybe,
  some,
  withDefault
} from './maybe';
import { components, dontDownloadForOutputTypes, DataEnvelope } from '../components/ui_components';
import { defaultApiHelpers } from '../contexts/api';

export function pick<T extends DataEnvelope<any>, K extends keyof T>(obj: T, keys: K[]) {
  return Object.fromEntries(
    keys.filter(key => key in obj).map(key => [key, obj[key]])
  );
};

export function pickToArray<T extends DataEnvelope<any>, K extends keyof T>(obj: T, keys: K[]) {
  return keys.filter(key => key in obj).map(key => obj[key]);
}

export function uniqueValues(original: any[]) {
  function onlyUnique(value: any, index: number, self: any) {
    return self.indexOf(value) === index;
  }
  return Array.from(original).filter(onlyUnique);
}

export function parseIfCan(data: any) {
  try {
    return JSON.parse(data);
  } catch {
    return data;
  }
}

export function unMaybe(obj: DataEnvelope<Maybe<any>>, keep: boolean = false): DataEnvelope<any> {
  // keeps only filled keys, and only their values
  if (keep) {
    // Not recommended due to CLASH: original null & [null] become the same. 
    return Object.fromEntries(
      Object.keys(obj).map((k) => [k, isSome(obj[k]) ? obj[k][0] : null])
    );
  }
  return Object.fromEntries(
    Object.keys(obj).filter((k) => isSome(obj[k])).map((k) => [k, obj[k][0]])
  );
}

function toSomes(obj: DataEnvelope<any>): DataEnvelope<Maybe<any>> {
  return Object.fromEntries(Object.keys(obj).map((k) => [k, some(obj[k])]))
}

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

// export function workflowById(
//   id: any,
//   state: VulcanState
// ): Workflow | undefined {
//   return state.workflows.find((w) => workflowId(w) === id);
// }

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

export function workspaceFromRaw(raw: WorkspaceRaw): Workspace {
  const vulcanConfig = vulcanConfigFromRaw(raw.vulcan_config);
  const intWorkspace: any = {...raw};
  delete intWorkspace['target_mapping']
  return {
    ...intWorkspace,
    steps: {...raw.target_mapping},
    vulcan_config: vulcanConfig
  };
}

export function compSetName(content: VulcanConfigElement) {
  // configUI
  if (!!content.output && !!content.output.params) {
    return content.output.params.join("/")
  }
  // inputUI && outputUI -- we'll require that at least one is given!
  return !!content.name ? content.name : content.display
}

// export function compSetLabel(content: VulcanConfigElement) {
//   return content.display
// }

export function vulcanConfigFromRaw(config: VulcanConfig | VulcanConfigRaw) {
  // Transforms it from array to object, using compSetNames as keys
  if (!Array.isArray(config)) return config;
  return Object.fromEntries(config.map(content => [compSetName(content), content]));
}

export function paramValuesFromRaw(param_values: Workspace['last_config'], workspace: Workspace): WorkspaceStatus['params'] {
  const {vulcan_config} = workspace;
  return Object.fromEntries(paramUINames(workspace).map(setName => [
      setName,
      !!param_values ?
        toSomes(pick(param_values, vulcan_config[setName].output?.params as string[])) :
        Object.fromEntries((vulcan_config[setName].output?.params as string[]).map(p => [p, null]))
    ]
  ))
}

export function paramValuesToRaw(params: WorkspaceStatus['params']): FlatParams {
  let output = {} as {[k: string]: any};
  Object.values(params).forEach(vals => {
    output = {...output, ...unMaybe(vals)};
  })
  return output;
}

export function uiContentsFromFiles(workspace: Workspace, file_contents?: WorkspaceStatus['file_contents']): WorkspaceStatus['ui_contents'] {
  const {vulcan_config} = workspace;
  if (!file_contents) {
    return Object.fromEntries(allUIStepNames(workspace).map(setName => [
      setName,
      Object.fromEntries((vulcan_config[setName].output?.files as string[]).map(f => [f, null]))
    ]));
  }
  return Object.fromEntries(allUIStepNames(workspace).map(setName => {
    const step_files = vulcan_config[setName].output?.files as string[];
    const fileData: DataEnvelope<Maybe<any>> = {};
    step_files.forEach(f => {
      if (f in file_contents) {
        fileData[f] = null;
      } else {
        fileData[f] = some(file_contents[f]);
      }
    })
    return [setName, fileData]
  }));
}

// // All components defined in the vulcan_config
// export function allUIComponentNames(
//   given: VulcanState | VulcanState['workspace']
// ): string[] {
//   const workspace = !!given && 'workspace' in given ? given.workspace : given;
//   return workspace ? Object.keys(workspace.vulcan_config) : [];
// }


// Components defined in the vulcan_config that populate "primary inputs" of the input_feed
//  - All will be ultimately lumped into the primary / initial input
//  - Uniquely: output to params / the config file
//  - Also WILL NOT exist in the 'workspace.dag' array of "steps"
export function paramUINames(
  given: VulcanState | VulcanState['workspace']
): string[] {
  const workspace = !!given && 'workspace' in given ? given.workspace : given;
  return !workspace ? [] : Object.keys(workspace.vulcan_config)
    .filter(k => !!workspace.vulcan_config[k].output && !!workspace.vulcan_config[k].output.params);
}

// Components in the vulcan_config that populate the input_feed and are also tracked as an explicit step of the workspace. 
//  - don't output as param / to the config.
//  - do output to a file.
//  - WILL exist in the 'workspace.dag' array of "steps" (generates an input for a later step, so tracked in the snakefile as a dummy step)
export function inputUINames(
  given: VulcanState | VulcanState['workspace']
): string[] {
  const workspace = !!given && 'workspace' in given ? given.workspace : given;
  if (!workspace) return [];
  let uiNames: string[] = [];
  for (const [key,value] of Object.entries(workspace.vulcan_config)) {
    // No params output, but has files output.  (MIGHT need to change this in the future)
    if (!!value.output && !value.output.params && value.output.files) uiNames.push(key);
  }
  return uiNames;
}

// Components in the vulcan_config that populate the output_feed
//  - don't have an output.
//  - do have an input from either param or file (more likely a file)
//  - SHOULD NOT exist in the 'workspace.dag' array of "steps", but could in the future
// MIGHT need to change this in the future for e.g. optional usage of a plot output by drawing a box to select observations for a followup.
export function outputUINames(
  given: VulcanState | VulcanState['workspace']
): string[] {
  const workspace = !!given && 'workspace' in given ? given.workspace : given;
  if (!workspace) return [];
  let uiNames: string[] = [];
  for (const [key,value] of Object.entries(workspace.vulcan_config)) {
    // No outputs.  (MIGHT need to change this in the future)
    if (!value.output) uiNames.push(key);
  }
  return uiNames;
}

export function allUIStepNames(
  given: VulcanState | VulcanState['workspace']
): string[] {
  const workspace = !!given && 'workspace' in given ? given.workspace : given;
  return workspace ? Object.keys(workspace.vulcan_config).filter(k => workspace.dag.includes(k)) : [];
}

export function uiNamesToBufferData(
  given: VulcanState | VulcanState['workspace']
): string[] {
  const workspace = !!given && 'workspace' in given ? given.workspace : given;
  if (!workspace) return [];
  // All inputUIs
  let uiNames: string[] = inputUINames(workspace);
  // Some outputUIs (based on the ui_component used)
  for (const out of outputUINames(workspace)) {
    if (!dontDownloadForOutputTypes.includes(workspace.vulcan_config[out].ui_component)) uiNames.push(out);
  }
  return uiNames;
}

export function allFilesToBuffer(
  given: VulcanState | VulcanState['workspace']
): string[] {
  const workspace = !!given && 'workspace' in given ? given.workspace : given;
  if (!workspace) return [];
  const {vulcan_config} = workspace;

  let inputFiles: string[] = [];
  for (let step in uiNamesToBufferData(workspace)) {
    inputFiles = inputFiles.concat(vulcan_config[step].output?.files as string[]);
  }
  return uniqueValues(inputFiles);
}

export function shouldDownloadOutput(
  workspace: Workspace,
  outputName: string
): boolean {
  const stepsToBufferFor = uiNamesToBufferData(workspace);
  const stepsAsInput = Object.keys(workspace.vulcan_config).filter(
    step => !!workspace.vulcan_config[step].output?.files && workspace.vulcan_config[step].output.files.includes(outputName)
  );
  return stepsAsInput.filter(step => stepsToBufferFor.includes(step)).length > 0;
}

export function filesReturnToMultiFileContent(filesContent: MultiFileContentResponse): MultiFileContent {
  const output: MultiFileContent = {};
  filesContent.forEach(f => {
    output[f.filename] = parseIfCan(f.content)
  });
  return output;
}

export function updateStepStatusesFromRunStatus(stepStatusReturns: RunStatus, stepStatusCurrent: WorkspaceStatus['steps']) {
  const newStepStatus = {...stepStatusCurrent};
  const newCompletions = [] as string[];
  Object.entries(stepStatusReturns).forEach(([stepName, statusFine]) => {
    if (newStepStatus[stepName].statusFine==statusFine) return
    const statusBroad = StatusStringBroaden(statusFine);
    newStepStatus[stepName]['status'] = statusBroad;
    newStepStatus[stepName]['statusFine'] = statusFine;
    if (statusBroad=='complete') {
      newCompletions.push(stepName);
    }
  });
  return {newStepStatus, newCompletions}
}

export function stepOfName(
  stepName: string,
  workspace_steps: Workspace['steps'],
  vulcan_config: VulcanConfig
): WorkspaceStep | undefined {
  if (stepName in workspace_steps) return workspace_steps[stepName];
  if (stepName in vulcan_config) {
    const step_conf = vulcan_config[stepName];
    return {
      name: stepName,
      input: step_conf.input || {},
      output: step_conf.output || {},
      vulcan_config: true,
      label: step_conf.display || stepName,
      ui_component: step_conf.ui_component,
      doc: step_conf.doc
    }
  }
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
  step: string | WorkspaceStep,
  status: VulcanState['status']
): StepStatus | undefined {
  const stepName = typeof step === 'string' ? step : step.name;
  return (stepName in status.steps) ? status.steps[stepName] : undefined;
}

export function statusStringOfStepOrGroupedStep(
  step: WorkspaceStep | WorkspaceStepGroup,
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
  step: WorkspaceStep | WorkspaceStepGroup
) {
  if ('steps' in step) {
    return step.label;
  }

  return step.label || step.name;
}

// ToDo: Assess need to keep
//  - additionally only outputs if an output component, but...
//  - could have a single component retrieval route for inputs and outputs
// export const uiOutputOfStep = (step: WorkspaceStep) => {
//   if (step.run && step.run.startsWith(RUN.UI_OUTPUT)) {
//     return step.run.split('/')[1].replace('.cwl', '');
//   }
//   return null;
// };

export const uiComponentOfStep = (step: string, vulcan_config: VulcanConfig) => {
  return (step in vulcan_config) ? vulcan_config[step].ui_component : null;
};

export const stepInputDataUrls = (
  step: WorkspaceStep,
  projectName: VulcanState['projectName'],
  workspaceId: VulcanState['workspaceId'],
  vulcanPath: typeof defaultApiHelpers.vulcanPath
): {[k: string]: string} => {
  // Parse input files of the step into a direct download url.
  const result: {[k: string]: string} = {};

  if (!!step.input.files) {
    step.input.files.forEach((file) => {
      result[file] = vulcanPath(`/api/v2/${projectName}/workspace/${workspaceId}/file/direct/${file}`)
    });
  }

  return result;
};

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

// export function compSetNameOfParam(paramName: string, vulcan_config: VulcanConfig) {
//   function hasParam(config: VulcanConfigElement) {
//     return !!config.output && !!config.output.params && config.output.params?.includes(paramName)
//   }
//   const found = Object.keys(vulcan_config).filter(
//     setName => hasParam(vulcan_config[setName])
//   )
//   return found.length > 0 ? found[0] : undefined;
// }

export function stepInputDataRaw(
  step: WorkspaceStep,
  last_params: WorkspaceStatus['last_params'],
  file_contents: WorkspaceStatus['file_contents']
): {[k: string]: any} {
  // Pull out any previous step's output data link that is a required
  //   input into this UI step.
  const result: {[k: string]: any} = {};

  if (!!step.input.params) {
    step.input.params.forEach((paramName) => {
      result[paramName] = last_params[paramName] || null;
    });
  }
  
  if (!!step.input.files) {
    step.input.files.forEach((fileName) => {
      result[fileName] = file_contents[fileName] || null;
    });
  }

  return result;
};

export function stepInputMapping(
  config: VulcanConfigElement,
) {
  // Returns the mapping of:
  //   - component's internal data value naming expectations ("internal-name")
  //   - actual file/param names to use as source ("external-name")
  // in the form:
  //   - {internal-name: external-name}
  const rawKeys = config.input_map ? config.input_map :
    config.input?.params ? config.input.params :
    config.input?.files ? config.input.files : [] as string[];
  // Todo MAYBE (should be enforced elsewhere): Error handling if empty array
  
  const componentNeeds = components[config.ui_component].slice(3,5) as [string[], string | undefined];
  // Required Elements
  const map: {[k: string]: string} = {};
  componentNeeds[0].forEach( (val, i) => {
    // ToDo MAYBE (should be enforced elsewhere): Error if rawKeys is too short
    map[val] = rawKeys[i]
  })

  // Optional Element
  if (rawKeys.length > componentNeeds[0].length && !!componentNeeds[1]) {
    map[componentNeeds[1]] = rawKeys[rawKeys.length-1];
  }

  return map;
}

export function stepOutputMapping(
  config: VulcanConfigElement,
) {
  // Returns the mapping of:
  //   - component's internal value naming expectations ("internal-name")
  //   - actual file/param names to output to ("external-name")
  // in the form:
  //   - {internal-name: external-name}
  const rawKeys = config.output?.params ? config.output.params :
    config.output?.files ? config.output.files : [] as string[];
  // Todo MAYBE (should be enforced elsewhere): Error handling if empty array
  
  const componentNeeds = components[config.ui_component][2];
  const map: {[k: string]: string} = {};
  componentNeeds.forEach( (val, i) => {
    // ToDo MAYBE (should be enforced elsewhere): Error if rawKeys is too short
    map[val] = rawKeys[i]
  })
  return map;
};

// export function stepOutputDataOriginal(
//   config: VulcanConfigElement,
//   last_params: WorkspaceStatus['last_params'],
//   file_contents: WorkspaceStatus['file_contents']
// ) {
//   const keyMap = stepOutputMapping(config);
//   const values: {[k: string]: any} = {};
//   Object.keys(keyMap).forEach( (name) => {
//     values[name] = name in last_params ?
//       last_params[keyMap[name]] :
//       name in file_contents ? file_contents[keyMap[name]] : null
//   });
//   return values;
// };

export function stepOutputData(
  stepName: string,
  keyMap: ReturnType<typeof stepOutputMapping>,
  buffered: DataEnvelope<Maybe<any>>,
  params: WorkspaceStatus['params'],
  ui_contents: WorkspaceStatus['ui_contents'],
  _default: any
) {
  const values: {[k: string]: Maybe<any>} = {};
  Object.entries(keyMap).forEach( ([internal, external], index) => {
    // internal = the value names known to the ui_component definition
    // external = the file or param name mapping to that value element
    const _def = index==0 ? _default : undefined;
    values[internal] = internal in buffered ? 
      buffered[internal] :
      stepName in params ?
      params[stepName][external] || _def || null :
      ui_contents[stepName][external] || _def || null
  });
  return values;
};

export function stepHasBeenOutput(
  stepName: string,
  status: WorkspaceStatus,
) {
  return (stepName in status.params) || (stepName in status.ui_contents);
}

export function pendingUIInputStepReady(
  step: string,
  status: VulcanState['status'],
  workspace: Workspace,
  data: VulcanState['data']
) {
  return (
    inputUINames(workspace).includes(step) &&
    pendingStepNames(workspace, status).includes(step) &&
    workspace.vulcan_config[step].input?.files?.every((id) => id in data)
  );
}

function stepNamesOfStatus(
  targetStatus: StatusString,
  workspace: Workspace,
  status: WorkspaceStatus
): string[] {
  return workspace.dag.filter(
    (step) => statusOfStep(step, status)?.status === targetStatus
  );
}
export function completedStepNames(workspace: Workspace, status: WorkspaceStatus): string[] {
  return stepNamesOfStatus('complete', workspace, status);
}
export function pendingStepNames(workspace: Workspace, status: VulcanState['status']): string[] {
  return stepNamesOfStatus('pending', workspace, status);
}
export function erroredStepNames(workspace: Workspace, status: VulcanState['status']): string[] {
  return stepNamesOfStatus('error', workspace, status);
}
export function upcomingStepNames(workspace: Workspace, status: VulcanState['status']): string[] {
  return stepNamesOfStatus('upcoming', workspace, status);
}

export function completedSteps(workspace: Workspace, status: WorkspaceStatus): WorkspaceStep[] {
  return pickToArray(workspace.steps, completedStepNames(workspace, status));
}
export function pendingSteps(workspace: Workspace, status: VulcanState['status']): WorkspaceStep[] {
  return pickToArray(workspace.steps, pendingStepNames(workspace, status));
}
export function erroredSteps(workspace: Workspace, status: VulcanState['status']): WorkspaceStep[] {
  return pickToArray(workspace.steps, erroredStepNames(workspace, status));
}

export function completedUiOutputSteps(
  workspace: Workspace,
  status: VulcanState['status']
): WorkspaceStep[] {
  const uiOutputSteps = outputUINames(workspace);
  const completedUIOutputStepNames = completedStepNames(workspace, status).filter(
    (step) => uiOutputSteps.includes(step)
  );
  return pickToArray(workspace.steps, completedUIOutputStepNames);
}

export function hasRunningSteps(status: VulcanState['status']): boolean {
  return Object.values(status.steps).filter((s) => s.status == 'running').length > 0;
}

export const inputGroupName = (name: string) => {
  let groupName = name.split('__')[0];
  if (groupName === name) return null;

  return groupName;
};

export function groupUiSteps(uiStepNames: string[], workspace: Workspace): WorkspaceStepGroup[] {
  const uiSteps = pickToArray(workspace.steps, uiStepNames);
  const map: {[k: string]: WorkspaceStepGroup} = {};
  const result: WorkspaceStepGroup[] = [];

  uiSteps.forEach((step) => {
    
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

export function useMemoized<P1, R>(f: (p1: P1) => R, p1: P1): R {
  return useMemo(() => f(p1), [f, p1]);
}
