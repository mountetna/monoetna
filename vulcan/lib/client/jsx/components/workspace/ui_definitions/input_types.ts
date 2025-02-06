import React, {useMemo} from 'react';
import {VulcanConfig, Workspace, WorkspaceStatus, WorkspaceStep} from '../../../api_types';
import {Maybe} from '../../../selectors/maybe';
import {
  paramUINames,
  stepInputDataRaw,
  stepInputMapping,
  stepOfName,
  stepOutputData,
  stepOutputMapping,
} from '../../../selectors/workflow_selectors';
import {defaultBufferedInputs} from '../../../contexts/input_state_management';
import {DataEnvelope} from '../../ui_components';
export {DataEnvelope} from '../../ui_components';

export function nulled_vals(de: DataEnvelope<any>): DataEnvelope<null> {
  let de_ = {...de};
  Object.keys(de_).forEach((i) => (de_[i] = null));
  return de_;
}

export interface InputSpecification {
  ui_component: string;
  label: string;
  name: string;
  doc?: string | null;
}

export interface BoundInputSpecification<Value = unknown, DataElement = unknown>
  extends InputSpecification {
  value: Maybe<Value>;

  onChange(v: Maybe<Value>, destructure?: boolean): void;

  data?: DataEnvelope<DataElement> | undefined | null;

  valueKeyMap: ReturnType<typeof stepOutputMapping>
}

export function getParamUISpecifications(
  workspace: Workspace | null
): InputSpecification[] {
  const elements = paramUINames(workspace);
  if (!workspace || elements.length < 1) return [];
  const {vulcan_config} = workspace;
  return elements.map(name => {
    const conf = {...vulcan_config[name]}
    return {
      ui_component: conf.ui_component,
      label: conf.display,
      name: name,
      doc: conf.doc || null
    }
  })
}

export function getInputSpecifications(
  step: WorkspaceStep | [string, WorkspaceStep][],
  workspace: Workspace | null
): InputSpecification[] {
  if (!workspace) return [];

  if (Array.isArray(step)) {
    return step.map(([name, WorkspaceStep]) => {
      return {
        ui_component: WorkspaceStep.ui_component || 'default',
        label: WorkspaceStep.label || name,
        name: name,
        doc: WorkspaceStep.doc || null
      };
    });
  }

  return [{
    name: step.name,
    doc: step.doc,
    label: step.label || step.name,
    ui_component: step.ui_component || 'default'
  }];
}

// export function collapseInputValues(
//   stepName: string | undefined,
//   inputName: string,
//   inputs:
//     | typeof defaultBufferedInputs.values
//     | VulcanState['session']['inputs'],
//   fromBuffer: boolean
// ): Maybe<unknown> {
//   const stepPrefixRegex = new RegExp(
//     stepName ? `^${inputName}\/` : `^${inputName}$`
//   );
//   const matchingInputs = Object.keys(inputs).filter((inputName: string) =>
//     inputName.match(stepPrefixRegex)
//   );

//   if (matchingInputs.length > 1) {
//     const groupedInputs = matchingInputs.reduce(
//       (acc: Maybe<{[key: string]: unknown}>, inputName) => {
//         if (!acc) return some({}); // should never return this? Here to make TSC happy.

//         const outputName = inputName.replace(stepPrefixRegex, '');
//         acc[outputName as any] = inputs[inputName];

//         return acc;
//       },
//       {} as Maybe<{[key: string]: unknown}>
//     );

//     return fromBuffer ? some(groupedInputs) : groupedInputs;
//   } else if (matchingInputs.length === 1) {
//     return inputs[matchingInputs[0]];
//   }

//   return null;
// }

export function bindInputSpecification(
  input: InputSpecification,
  // workspace_steps: Workspace['steps'],
  vulcan_config: VulcanConfig,
  last_params: WorkspaceStatus['last_params'],
  file_contents: WorkspaceStatus['file_contents'],
  params: WorkspaceStatus['params'],
  ui_contents: WorkspaceStatus['ui_contents'],
  buffered: typeof defaultBufferedInputs.values,
  setValues: typeof defaultBufferedInputs.setValues
): BoundInputSpecification {
  
  const stepName = input.name;
  const step = stepOfName(stepName, 
    // workspace_steps,
    vulcan_config);
  // ToDo: Surface an error instead
  if (!step) return;
  const config = vulcan_config[stepName];

  const inputDataKeyMap = stepInputMapping(config);
  const inputDataRaw = stepInputDataRaw(step, last_params, file_contents);
  const inputData: {[k: string]: any} = {};
  Object.keys(inputDataKeyMap).forEach( (name) => {
    inputData[name] = inputDataRaw[inputDataKeyMap[name]];
  })

  const outputDataKeyMap = stepOutputMapping(config);
  const value = stepOutputData(stepName, outputDataKeyMap, buffered, params, ui_contents, config['default']);

  return {
    ...input,
    data: inputData,
    value: value,
    onChange(v: {[k:string]: Maybe<unknown>}) {
      setValues(v);
    },
    valueKeyMap: outputDataKeyMap
  };
}

export type WorkspaceStepGroup = {label: string; steps: WorkspaceStep[]};

export type InputBackendComponent<
  Params extends {} = {},
  Value = unknown,
  DataElement = unknown
> = (
  p: WithInputParams<Params, Value, DataElement>
) => React.ReactElement | null;
export type WithInputParams<
  Params extends {},
  Value,
  DataElement = unknown
> = Params & {
  onChange: (v: Maybe<Value>, destructure?: boolean) => void;
  value: Maybe<Value>;
  data: DataEnvelope<DataElement> | undefined | null;
  numOutputs?: number;
};

export interface ValidationInputSpecification<
  Value = unknown,
  DataElement = unknown
> {
  data?: DataEnvelope<DataElement> | null;
  value: Maybe<Value>;
}

export type InputValidator<Value = unknown, DataElement = unknown> = (
  input: ValidationInputSpecification<Value, DataElement>
) => string[];
