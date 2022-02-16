import * as React from 'react';
import {Workflow, WorkflowInput, WorkflowStep} from '../../../../api_types';
import {Maybe, some} from "../../../../selectors/maybe";
import {Dispatch} from "react";
import {VulcanAction} from "../../../../actions/vulcan_actions";
import {VulcanState} from "../../../../reducers/vulcan_reducer";
import {
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
  console.log('buffered', buffered, input, session.inputs)
  return {
    ...input,
    onChange(v: Maybe<unknown>) {
      setInputs(inputs => ({...inputs, [input.source]: v}))
    },
    data: inputData,
    value: input.source in buffered ? buffered[input.source] :
      input.source in session.inputs ? some(session.inputs[input.source]) :
        null,
  };
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
