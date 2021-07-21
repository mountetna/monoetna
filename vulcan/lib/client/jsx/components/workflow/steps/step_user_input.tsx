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
import {useWorkflow} from "../../../contexts/workflow_context";

export default function StepUserInput({ step, }: { step: WorkflowStep; }) {
  const {state, dispatch} = useContext(VulcanContext);
  const {workflow} = useWorkflow();
  const uiQuery = uiQueryOfStep(step);

  const stepInputs: BoundInputSpecification[] = useMemo(() => getInputSpecifications(step, workflow)
    .map(spec => bindInputSpecification(spec, workflow, state.status, state.session, state.data, state.bufferedInputValues, dispatch)),
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
