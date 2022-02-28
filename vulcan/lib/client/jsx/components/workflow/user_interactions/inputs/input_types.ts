import * as React from 'react';
import {Workflow, WorkflowInput, WorkflowStep} from '../../../../api_types';
import {Maybe, some, isSome, withDefault} from "../../../../selectors/maybe";
import {Dispatch} from "react";
import {VulcanAction} from "../../../../actions/vulcan_actions";
import {VulcanState} from "../../../../reducers/vulcan_reducer";
import {
  sourceNamesOfStep, splitSource, stepInputDataRaw, stepOfSource, stepOfStatus, stepOutputs, uiQueryOfStep
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

  onChange(v: Maybe<Value>, destructure?: boolean): void;

  data?: DataEnvelope<DataElement> | undefined | null;

  numOutputs?: number;
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
    onChange(v: Maybe<unknown>, destructure: boolean = false) {
      if (destructure) {
        const authoredOutputs = stepOutputs(step);
        
        let values: {[key: string]: any} = {};
        const userValue: {[key: string]: any} = withDefault(v, {}) as any;
        // First, look for authored outputs that appear in the userValue object
        const [withValues, withoutValues] = authoredOutputs.reduce((acc, authoredOutput) => {
          let [_, outputName] = splitSource(authoredOutput);
          
          acc[userValue.hasOwnProperty(outputName) ? 0 : 1].push(authoredOutput);

          return acc;
        }, [[] as string[], [] as string[]])
        
        withValues.forEach((authoredOutput) => {
          let [_, outputName] = splitSource(authoredOutput);
            // Widget has specified a Hash that explicitly matches
            //   the CWL output names.
            values[authoredOutput] = userValue[outputName as any];
        })

        // Next, any remaining, unassigned outputs are picked
        //   from the unassigned userValue keys as if we were popping
        //   values off a queue.
        const usedValueKeys = Object.keys(values).map((k) => splitSource(k)[1]);
        const unassignedValueKeys = Object.keys(userValue).filter((key: string) => {
          return !usedValueKeys.includes(key)
        }).sort();
        
        withoutValues.forEach((authoredOutput, index) => {
          if (unassignedValueKeys[index])
            values[authoredOutput] = userValue[unassignedValueKeys[index] as any];
        });

        // Any unassigned values in `userValue` are discarded at this point.
        // Warn in the console. Should we warn the user more visibly?
        if (unassignedValueKeys.length > withoutValues.length) {
          console.warn("UI input returned values not assigned in CWL.", unassignedValueKeys.slice(withoutValues.length))
        }
        
        setInputs(inputs => ({...inputs, ...values as any}))
      } else {
        setInputs(inputs => ({...inputs, [input.source]: v}))
      }
    },
    data: inputData,
    value: input.source in buffered ?
          collapseInputValues(input.name, buffered, true) :
      input.source in session.inputs ?
        some(collapseInputValues(input.name, session.inputs, false)) :
      null,
    numOutputs: stepOutputs(stepOfStatus(input.name, workflow)).length
  };
}

function collapseInputValues(
  stepName: string,
  inputs: typeof defaultBufferedInputs.inputs | VulcanState['session']['inputs'],
  fromBuffer: boolean
): Maybe<unknown> {
  const stepPrefixRegex = new RegExp(`^${stepName}\/?`);
  const matchingInputs = Object.keys(inputs).filter((inputName: string) =>
    inputName.match(stepPrefixRegex)
  );

  if (matchingInputs.length > 1) {
    const groupedInputs = matchingInputs.reduce((acc: Maybe<{[key: string]: unknown}>, inputName) => {
      if (!acc) return some({}); // should never return this? Here to make TSC happy.

      const outputName = inputName.replace(stepPrefixRegex, '');
      acc[outputName as any] = inputs[inputName];

      return acc;
    }, {} as Maybe<{[key: string]: unknown}>);
    return fromBuffer ? some(groupedInputs) : groupedInputs;
  } else if (matchingInputs.length === 1) {
    return inputs[matchingInputs[0]];
  }

  return null;
}

export type WorkflowStepGroup = { label: string, steps: WorkflowStep[] }

export type InputBackendComponent<Params extends {} = {}, Value = unknown, DataElement = unknown> = (p: WithInputParams<Params, Value, DataElement>) => React.ReactElement | null;
export type WithInputParams<Params extends {}, Value, DataElement = unknown> = Params & {
  onChange: (v: Maybe<Value>, destructure?: boolean) => void;
  value: Maybe<Value>;
  data: DataEnvelope<DataElement> | undefined | null;
  numOutputs?: number;
}

export interface ValidationInputSpecification<Value = unknown, DataElement = unknown> {
  data?: DataEnvelope<DataElement> | null;
  value: Maybe<Value>,
}

export type InputValidator<Value = unknown, DataElement = unknown> = (input: ValidationInputSpecification<Value, DataElement>) => string[];
