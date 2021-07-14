import React, { useContext, useMemo } from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import UserInput from '../user_interactions/inputs/user_input';

import {WorkflowStep} from '../../../api_types';
import {
  uiQueryOfStep,
} from '../../../selectors/workflow_selectors';
import {
  bindInputSpecification, BoundInputSpecification, getInputSpecifications
} from '../user_interactions/inputs/input_types';

export default function StepUserInput({ step, }: { step: WorkflowStep; }) {
  const {state, dispatch} = useContext(VulcanContext);
  const {workflow} = state;
  const uiQuery = uiQueryOfStep(step);

  const stepInputs: BoundInputSpecification[] = useMemo(() => getInputSpecifications(step, workflow)
    .map(spec => bindInputSpecification(spec, state, dispatch)),
    [step, workflow, state, dispatch])

  if (!uiQuery) return null;

  return (
    <React.Fragment>
      {stepInputs.map((input, index) => {
        return (
          <UserInput
            input={input}
            hideLabel={true}
            key={index}
          />
        );
      })}
    </React.Fragment>
  );
}
