import React, { useContext, useMemo } from 'react';
import {VulcanContext} from '../../../contexts/vulcan_context';
import UserInput from '../user_interactions/inputs/user_input';
import {WorkflowStep} from '../../../api_types';
import {
  bindInputSpecification, BoundInputSpecification, getInputSpecifications, InputSpecification
} from '../user_interactions/inputs/input_types';
import {useWorkflow} from "../../../contexts/workflow_context";
import {BufferedInputsContext, WithBufferedInputs} from "../../../contexts/input_state_management";

export default function StepUserInput({ step, hideLabel = true }: { step: WorkflowStep; hideLabel: boolean }) {
  const {dispatch, commitSessionInputChanges} = useContext(VulcanContext);
  const {workflow} = useWorkflow();
  const specs = useMemo(() => getInputSpecifications(step, workflow), [step, workflow]);

  return (
    <React.Fragment>
      <WithBufferedInputs commitSessionInputChanges={commitSessionInputChanges} dispatch={dispatch} stepName={step.name}>
        <StepUserInputInner specs={specs} hideLabel={hideLabel}/>
      </WithBufferedInputs>
    </React.Fragment>
  );
}

function StepUserInputInner({ specs, hideLabel }: {specs: InputSpecification[], hideLabel: boolean }) {
  const {state} = useContext(VulcanContext);
  const {status, session, data} = state;
  const {workflow} = useWorkflow();
  const {inputs, setInputs} = useContext(BufferedInputsContext);

  const stepInputs: BoundInputSpecification[] = useMemo(() => specs
    .map(spec => bindInputSpecification(spec, workflow, status, session, data, inputs, setInputs)),
    [specs, workflow, status, session, data, inputs, setInputs])

  return (
    <React.Fragment>
      {stepInputs.map((input, index) => {
        return (
          <UserInput
            input={input}
            hideLabel={hideLabel}
            key={index}
          />
        );
      })}
    </React.Fragment>
  );
}
