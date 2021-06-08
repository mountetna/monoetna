import React, {
  useContext,
  useState,
  useEffect,
  useCallback,
  useMemo
} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import UserInput from '../user_interactions/inputs/user_input';

import {WorkflowStep} from '../../../api_types';
import {
  sourceNamesOfStep,
  uiQueryOfStep,
  stepInputDataRaw
} from '../../../selectors/workflow_selectors';
import {
  InputSpecification,
  InputType
} from '../user_interactions/inputs/input_types';

export default function StepUserInput({
  step,
  handleInputChange
}: {
  step: WorkflowStep;
  handleInputChange: (source: string, val: any) => void;
}) {
  const {state} = useContext(VulcanContext);
  const uiQuery = uiQueryOfStep(step);

  // We need to map the user input step's output to
  //   a set of input items.
  const outputRefs = useMemo(() => sourceNamesOfStep(step), [step]);
  const stepInputs: InputSpecification[] = useMemo(() => outputRefs.map(outputName => ({
    type: uiQuery as InputType,
    label: step.label || step.name,
    // The existing value
    value: state.inputs[outputName],
    name: outputName,
    data: stepInputDataRaw(step, state.status, state.data, state.session),
    doc: step.doc
  })), [outputRefs, uiQuery, step, state.inputs, state.status, state.data, state.session])

  if (!uiQuery) return null;

  return (
    <React.Fragment>
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
    </React.Fragment>
  );
}
