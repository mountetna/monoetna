import React, {useContext, useMemo} from 'react';
import {VulcanContext} from '../../../contexts/vulcan_context';
import UserInput from '../user_interactions/inputs/user_input';
import {
  GroupedInputStep, BoundInputSpecification, getInputSpecifications, bindInputSpecification
} from '../user_interactions/inputs/input_types';
import {sortInputsByLabel} from '../../../selectors/workflow_selectors';
import {useWorkflow} from "../../../contexts/workflow_context";

export default function StepUserInputDrawer({step}: { step: GroupedInputStep; }) {
  let {state, dispatch} = useContext(VulcanContext);
  const {workflow} = useWorkflow();

  // We need to unpack the grouped steps and add docs
  let stepInputs: BoundInputSpecification[] = useMemo(() => getInputSpecifications(step, workflow)
    .map(spec => bindInputSpecification(spec,
      workflow,
      state.status,
      state.session,
      state.data,
      state.bufferedInputValues,
      dispatch
    )), [step, workflow, state, dispatch]);

  if (!workflow) return null;

  return (<React.Fragment>
    {sortInputsByLabel(stepInputs).map((input, index) => {
      return (<UserInput input={input} hideLabel={false} key={index}/>);
    })}
  </React.Fragment>);
}
