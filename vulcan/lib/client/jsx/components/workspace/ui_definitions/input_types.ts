import React, {useCallback, useMemo} from 'react';
import {VulcanConfig, VulcanConfigElement, InputConfig, Workspace, WorkspaceStatus, WorkspaceStep} from '../../../api_types';
import {Maybe} from '../../../selectors/maybe';
import {
  paramUINames,
  stepOfName,
} from '../../../selectors/workflow_selectors';
import {defaultBufferedInputs} from '../../../contexts/input_state_management';
import {DataEnvelope, inputComponents} from '../../ui_components';

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

  onChange(v: Maybe<Value>): void;

  data?: DataEnvelope<DataElement> | undefined | null;

  valueKeyMap: ReturnType<typeof stepOutputMapping>;

  defaultValue?: any;

  showError: (e: string) => void;
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

export function fillInputData(
  config: VulcanConfigElement,
  step: WorkspaceStep,
  last_params: WorkspaceStatus['last_params'],
  file_contents: WorkspaceStatus['file_contents'],
  showError?: ((e: string) => void)
) {

  // First: what we have from the vulcan_config
  const mappedConfigInputs = config.input?.preset || {}
  let unmappedConfigInputs: any[] = [];

  if (!!step.input.params && !!step.input.files && !!showError) {
    showError("Workflow Issue: VulcanConfigElements can have params or files inputs, but NOT both.")
  }
  function fill_input_key_from_data(
    def: InputConfig['params'] | InputConfig['files'],
    currentVals: WorkspaceStatus['last_params'] | WorkspaceStatus['file_contents']
  ) {
    if (!!def) {
      if (Array.isArray(def)) {
        def.forEach((name) => {
          unmappedConfigInputs.push(currentVals[name]);
        });
      } else {
        Object.keys(def).forEach((name) => {
          if (name in mappedConfigInputs && !!showError) {
            showError("Workflow Issue: duplicate input names in VulcanConfigElement")
          };
          mappedConfigInputs[name] = currentVals[def[name]]
        });
      }
    }
  }
  fill_input_key_from_data(step.input.params, last_params)
  fill_input_key_from_data(step.input.files, file_contents)

  // Now: build with context from the ui_component
  if (!(config.ui_component in inputComponents)) {
    // stepOutput UI
    unmappedConfigInputs.forEach((value) => mappedConfigInputs[value] = value);
    return mappedConfigInputs;
  }
  // stepInput or Param UI
  const componentNeeds = inputComponents[config.ui_component][3] as string[];
  const componentOptional = inputComponents[config.ui_component][4] as string | undefined;

  const unusedKeys = Object.keys(mappedConfigInputs).filter((name) => !componentNeeds.includes(name) && name != componentOptional);
  if (unusedKeys.length > 0 && !!showError) {
    showError(`Specified input key(s) '${unusedKeys.join("', '")}' in VulcanConfigElement do not exist for ui_component '${config.ui_component}'`)
  }

  // Required Elements
  const data: {[k: string]: any} = {};
  componentNeeds.forEach( (comp) => {
    // ToDo: Error if nothing unmapped left to grab
    data[comp] = comp in mappedConfigInputs ?
      mappedConfigInputs[comp] :
      unmappedConfigInputs.shift()
  })

  // Optional Element
  if (!!componentOptional) {
    if (unmappedConfigInputs.length > 0 || componentOptional in mappedConfigInputs) {
      data[componentOptional] = componentOptional in mappedConfigInputs ?
        mappedConfigInputs[componentOptional] :
        unmappedConfigInputs.shift();
    }
  }

  return data
};

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
  
  const componentNeeds = inputComponents[config.ui_component][2];
  const map: {[k: string]: string} = {};
  componentNeeds.forEach( (val, i) => {
    // ToDo MAYBE (should be enforced elsewhere): Error if rawKeys is too short
    map[val] = rawKeys[i]
  })
  return map;
};

function stepOutputData(
  stepName: string,
  keyMap: ReturnType<typeof stepOutputMapping>,
  buffered: DataEnvelope<Maybe<any>>,
  params: WorkspaceStatus['params'],
  ui_contents: WorkspaceStatus['ui_contents']
) {
  const values: {[k: string]: Maybe<any>} = {};
  Object.entries(keyMap).forEach( ([internal, external], index) => {
    // internal = the value names known to the ui_component definition
    // external = the file or param name mapping to that value element
    values[internal] = internal in buffered ? 
      buffered[internal] :
      stepName in params ?
      params[stepName][external] || null :
      ui_contents[stepName][external] || null
  });
  return values;
};

export function bindInputSpecification(
  input: InputSpecification,
  // workspace_steps: Workspace['steps'],
  vulcan_config: VulcanConfig,
  last_params: WorkspaceStatus['last_params'],
  file_contents: WorkspaceStatus['file_contents'],
  params: WorkspaceStatus['params'],
  ui_contents: WorkspaceStatus['ui_contents'],
  buffered: typeof defaultBufferedInputs.values,
  setValues: typeof defaultBufferedInputs.setValues,
  showError: typeof defaultBufferedInputs.showError
): BoundInputSpecification {
  
  const stepName = input.name;
  const step = stepOfName(stepName, vulcan_config);
  // ToDo: Surface an error instead
  const config = vulcan_config[stepName];

  const outputDataKeyMap = stepOutputMapping(config);
  const value = stepOutputData(stepName, outputDataKeyMap, buffered, params, ui_contents);

  return {
    ...input,
    data: fillInputData(config, step, last_params, file_contents, showError),
    value: value,
    onChange(v: {[k:string]: Maybe<unknown>}) {
      setValues(v);
    },
    valueKeyMap: outputDataKeyMap,
    defaultValue: config['default'],
    showError
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
  defaultValue?: any;
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
