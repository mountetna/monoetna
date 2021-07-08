import * as React from 'react';
import {WorkflowStep} from '../../../../api_types';
import {Maybe} from "../../../../selectors/maybe";
import {Dispatch} from "react";

export type InputType = string;

export interface InputSpecification {
  type: InputType;
  label: string;
  value: Maybe<any>;
  onChange(v: Maybe<any>): void;
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

export type InputBackendComponent = (p: WithInputParams<{}>) => React.ReactElement | null;
export type WithInputParams<T extends {}> = T & {
  onChange: InputSpecification['onChange'];
  value: InputSpecification['value'];
  data: InputSpecification['data'];
}

export type InputValidator = (input: InputSpecification) => string[];
