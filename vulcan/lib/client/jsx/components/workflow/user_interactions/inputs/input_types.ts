import * as React from 'react';
import {WorkflowStep} from '../../../../api_types';
import {Maybe} from "../../../../selectors/maybe";

export type InputType = string;

// ui-queries using CWL take a dictionary of inputs, which often must be flattened.
// This type represents the unflattened raw data format coming in from the cwl.
export type DataEnvelope<Inner> = {[k: string]: Inner};

export interface InputSpecification<Value = unknown, DataElement = unknown> {
  type: InputType;
  label: string;
  value: Maybe<Value>;
  onChange(v: Maybe<Value>): void;
  data?: DataEnvelope<DataElement> | undefined | null;
  doc?: string | null;
}

export interface UIStep {
  step: WorkflowStep | GroupedInputStep;
  // Index in the original workflow
  index: number;
}

export type GroupedInputStep = WorkflowStep & {
  isGroup: true;
};

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
