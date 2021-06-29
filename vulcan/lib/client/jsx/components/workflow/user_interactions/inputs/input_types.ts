import * as React from 'react';
import {WorkflowStep} from '../../../../api_types';

export type InputType = string;

export interface InputSpecification {
  type: InputType;
  label: string;
  value: null | [any];
  data?: {[k: string]: any} | null;
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

export type InputOnChange = (inputName: string, val: any) => void;

export type InputBackendComponent = (p: {
  input: InputSpecification;
  onChange: InputOnChange;
  onClear?: () => void;
  onAll?: () => void;
}) => React.ReactElement | null;

export type InputValidator = (input: InputSpecification) => string[];
