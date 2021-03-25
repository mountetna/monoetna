import React, {useContext, useState, useEffect, useCallback, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {STATUS} from '../../../api_types';
import UserInput from '../user_interactions/inputs/user_input';
import StepName from './step_name';

import {WorkflowStep} from "../../../api_types";
import {
  sourceNamesOfStep,
  statusOfStep,
  uiQueryOfStep, stepInputDataRaw,
} from "../../../selectors/workflow_selectors";
import {InputSpecification, InputType} from "../user_interactions/inputs/input_types";

export default function StepUserInput({step, handleInputChange}: {step: WorkflowStep, handleInputChange: (source: string, val: any) => void}) {
  const [open, setOpen] = useState(true);
  const {state} = useContext(VulcanContext);

  const status = statusOfStep(step, state.status);
  const uiQuery = uiQueryOfStep(step);

  if (!status || !uiQuery) return null;

  const toggleInputs = useCallback(() => setOpen(!open), [setOpen, open]);

  // We need to map the user input step's output to
  //   a set of input items.
  const outputRefs = useMemo(() => sourceNamesOfStep(step), [step]);
  const stepInputs: InputSpecification[] = useMemo(() => outputRefs.map(outputName => ({
    type: uiQuery as InputType,
    label: step.label || step.name,
    // The existing value
    default: state.inputs[outputName],
    name: outputName,
    data: stepInputDataRaw(step, state.status, state.data, state.session),
    doc: step.doc
  })), [outputRefs, step, state.status, state.data, state.inputs])

  useEffect(() => {
    let stepStatus = status.status;
    if (STATUS.COMPLETE === stepStatus) {
      if (stepInputs.every(({data}) => !!data))
        setOpen(false);
    } else if (STATUS.PENDING === stepStatus) {
      setOpen(true);
    }
  }, [status.status]);

  return (
    <div className='step-user-input step'>
      <div onClick={toggleInputs}>
        <StepName step={step} />
      </div>
      <div
        className={`step-user-input-inputs sliding-panel vertical ${
          open ? 'open' : 'closed'
        }`}
      >
        {stepInputs.map((input, index) => {
          return (
            <UserInput
              input={input}
              hideLabel={true}
              onChange={handleInputChange}
              key={index}
            />
          );
        })}
      </div>
    </div>
  );
}
