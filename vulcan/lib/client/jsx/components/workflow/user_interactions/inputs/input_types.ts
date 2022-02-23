import * as React from 'react';
import {Workflow, WorkflowInput, WorkflowStep} from '../../../../api_types';
import {Maybe, some, isSome, withDefault} from "../../../../selectors/maybe";
import {Dispatch} from "react";
import {VulcanAction} from "../../../../actions/vulcan_actions";
import {VulcanState} from "../../../../reducers/vulcan_reducer";
import {
  collapsesOutputs,
  sourceNamesOfStep, stepInputDataRaw, stepOfSource, stepOfStatus, uiQueryOfStep
} from "../../../../selectors/workflow_selectors";
import {defaultBufferedInputs} from "../../../../contexts/input_state_management";

export type InputType = string;

// ui-queries using CWL take a dictionary of inputs, which often must be flattened.
// This type represents the unflattened raw data format coming in from the cwl.
export type DataEnvelope<Inner> = { [k: string]: Inner };

export function nulled_vals(de: DataEnvelope<any>): DataEnvelope<null> {
  let de_ = {...de};
  Object.keys(de_).forEach((i) => de_[i] = null);
  return de_
}

export interface InputSpecification {
  type: InputType;
  label: string;
  name: string;
  doc?: string | null;
  source: string,
}

export interface BoundInputSpecification<Value = unknown, DataElement = unknown> extends InputSpecification {
  value: Maybe<Value>;

  onChange(v: Maybe<Value>): void;

  data?: DataEnvelope<DataElement> | undefined | null;
}

export function getInputSpecifications(
  step: WorkflowStep | [string, WorkflowInput][],
  workflow: Workflow | null,
): InputSpecification[] {
  if (!workflow) return [];

  if (Array.isArray(step)) {
    return step.map(([name, workflowInput]) => {
      return {
        ...workflowInput,
        name,
        label: workflowInput.label || name,
        source: name,
      }
    })
  }

  return sourceNamesOfStep(step).map(source => ({
    source, name: step.name, doc: step.doc, label: step.label || step.name, type: uiQueryOfStep(step) || 'default'
  }));
}

export function bindInputSpecification(input: InputSpecification,
  workflow: Workflow,
  status: VulcanState['status'],
  session: VulcanState['session'],
  data: VulcanState['data'],
  buffered: typeof defaultBufferedInputs.inputs,
  setInputs: typeof defaultBufferedInputs.setInputs,
): BoundInputSpecification {
  const stepName = stepOfSource(input.source);
  const step = stepName && stepOfStatus(stepName, workflow);
  const inputData = stepName && step && stepInputDataRaw(step, status, data, session) || {};

  return {
    ...input,
    onChange(v: Maybe<unknown>) {
      if (collapsesOutputs(input.type) && isSome(v)) {
        const values = Object.entries(withDefault(v, {}) as any).reduce((acc: {[key: string]: unknown}, [key, value]) => {
          acc[`${input.name}/${key}`] = value;

          return acc;
        }, {});
        setInputs(inputs => ({...inputs, ...values as any}))
      } else {
        setInputs(inputs => ({...inputs, [input.source]: v}))
      }
    },
    data: inputData,
    value: input.source in buffered ?
        collapsesOutputs(input.type) ?
          collapseInputValues(input.name, buffered) :
          buffered[input.source] :
      input.source in session.inputs ?
        collapsesOutputs(input.type) ?
          collapseInputValues(input.name, session.inputs) :
          some(session.inputs[input.source]) :
        null,
  };
}

function collapseInputValues(stepName: string, inputs: typeof defaultBufferedInputs.inputs | VulcanState['session']['inputs']): Maybe<unknown> {
  const stepPrefix = `${stepName}/`;
  return some(Object.keys(inputs).filter((inputName: string) => {
    return inputName.startsWith(stepPrefix)
  }).reduce((acc: {[key: string]: unknown}, inputName) => {
    const outputName = inputName.replace(stepPrefix, '');
    acc[outputName as any] = inputs[inputName];
    return acc;
  }, {}));
}

export type WorkflowStepGroup = { label: string, steps: WorkflowStep[] }

export type InputBackendComponent<Params extends {} = {}, Value = unknown, DataElement = unknown> = (p: WithInputParams<Params, Value, DataElement>) => React.ReactElement | null;
export type WithInputParams<Params extends {}, Value, DataElement = unknown> = Params & {
  onChange: (v: Maybe<Value>) => void;
  value: Maybe<Value>;
  data: DataEnvelope<DataElement> | undefined | null;
}

export interface ValidationInputSpecification<Value = unknown, DataElement = unknown> {
  data?: DataEnvelope<DataElement> | null;
  value: Maybe<Value>,
}

export type InputValidator<Value = unknown, DataElement = unknown> = (input: ValidationInputSpecification<Value, DataElement>) => string[];
