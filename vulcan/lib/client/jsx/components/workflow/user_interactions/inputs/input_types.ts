import * as React from 'react';
import {defaultWorkflow, defaultWorkflowStep, Workflow, WorkflowInput, WorkflowStep} from '../../../../api_types';
import {Maybe} from "../../../../selectors/maybe";
import {Dispatch} from "react";
import {setBufferedInput, VulcanAction} from "../../../../actions/vulcan_actions";
import {VulcanState} from "../../../../reducers/vulcan_reducer";
import {
  sourceNamesOfStep, stepInputDataRaw, stepOfSource, stepOfStatus, uiQueryOfStep
} from "../../../../selectors/workflow_selectors";

export type InputType = string;

// ui-queries using CWL take a dictionary of inputs, which often must be flattened.
// This type represents the unflattened raw data format coming in from the cwl.
export type DataEnvelope<Inner> = { [k: string]: Inner };

export interface InputSpecification {
  type: InputType;
  label: string;
  doc?: string | null;
  source: string,
}

export interface BoundInputSpecification<Value = unknown, DataElement = unknown> extends InputSpecification {
  value: Maybe<Value>;

  onChange(v: Maybe<Value>): void;

  data?: DataEnvelope<DataElement> | undefined | null;
}

export function getInputSpecifications(
  step: WorkflowStep | GroupedInputStep | [string, WorkflowInput][],
  workflow: Workflow | null,
): InputSpecification[] {
  if (!workflow) return [];

  if ('isGroup' in step) {
    return step.in.map(({source, label, doc}) => ({
      source,
      doc,
      label: label || source,
      type: uiQueryOfStep(stepOfStatus(stepOfSource(source) || '', workflow) || defaultWorkflowStep) || 'default'
    }))
  }

  if (Array.isArray(step)) {
    return step.map(([name, workflowInput]) => {
      return {
        ...workflowInput,
        label: workflowInput.label || name,
        source: name,
      }
    })
  }

  return sourceNamesOfStep(step).map(source => ({
    source, doc: step.doc, label: step.label || step.name, type: uiQueryOfStep(step) || 'default'
  }));
}

export function bindInputSpecification(input: InputSpecification,
  workflow: Workflow,
  status: VulcanState['status'],
  session: VulcanState['session'],
  data: VulcanState['data'],
  buffered: VulcanState['bufferedInputValues'],
  dispatch: Dispatch<VulcanAction>
): BoundInputSpecification {
  const stepName = stepOfSource(input.source);
  const step = stepName && stepOfStatus(stepName, workflow);
  const inputData = stepName && step && stepInputDataRaw(step, status, data, session) || {};

  return {
    ...input,
    onChange(v: Maybe<unknown>) {
      dispatch(setBufferedInput({...buffered, [input.source]: v}))
    },
    data: inputData,
    value: input.source in buffered ?
      buffered[input.source] :
      input.source in session.inputs ? session.inputs[input.source] : null,
  };
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
  onChange: (v: Maybe<Value>) => void; value: Maybe<Value>; data: DataEnvelope<DataElement> | undefined | null;
}

export interface ValidationInputSpecification<Value = unknown, DataElement = unknown> {
  data?: DataEnvelope<DataElement> | null;
  value: Maybe<Value>,
}

export type InputValidator<Value = unknown, DataElement = unknown> = (input: ValidationInputSpecification<Value, DataElement>) => string[];
