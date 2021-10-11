import React, {useMemo} from 'react';
import { WorkflowStepGroup } from '../user_interactions/inputs/input_types';
import StepUserInput from "./step_user_input";

const collator = new Intl.Collator(undefined, {
  numeric: true, sensitivity: 'base'
});

export default function StepUserInputDrawer({group}: { group: WorkflowStepGroup; }) {
  const {steps} = group;

  // We need to unpack the grouped steps and add docs
  let stepInputs = useMemo(() => steps.sort(
    (a, b) => collator.compare(a.label || a.name, b.label || b.name))
    .map(step => <StepUserInput key={step.name} step={step} hideLabel={true}/>), [steps]);

  return (<React.Fragment>
    {stepInputs}
  </React.Fragment>);
}
