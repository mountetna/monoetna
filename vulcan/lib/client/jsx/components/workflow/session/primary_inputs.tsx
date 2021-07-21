import React, {useCallback, useContext, useEffect, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {inputGroupName} from '../../../selectors/workflow_selectors';
import InputGroup from './input_group';
import {
  bindInputSpecification,
  BoundInputSpecification, getInputSpecifications
} from "../user_interactions/inputs/input_types";
import {useWorkflow} from "../../../contexts/workflow_context";
import {maybeOfNullable} from "../../../selectors/maybe";
import {setBufferedInput} from "../../../actions/vulcan_actions";

export default function PrimaryInputs() {
  const {state, dispatch} = useContext(VulcanContext);
  const {bufferedInputValues, session: {inputs}} = state;
  const {workflow} = useWorkflow();

  // Ensure defaults are set.
  useEffect(() => {
    let updatedInputs = bufferedInputValues;
    Object.keys(workflow.inputs).forEach(inputName => {
      if (!(inputName in inputs) && !(inputName in bufferedInputValues)) {
        if (updatedInputs === bufferedInputValues) updatedInputs = {...updatedInputs};
        updatedInputs[inputName] = maybeOfNullable(workflow.inputs[inputName].default);
      }
    })

    if (updatedInputs !== bufferedInputValues) {
      dispatch(setBufferedInput(updatedInputs));
    }
  }, [bufferedInputValues, dispatch, inputs, workflow.inputs])

  const inputSpecifications = useMemo(() =>
    getInputSpecifications(Object.entries(workflow.inputs), workflow), [workflow]);

  let groupedInputs = useMemo(() => {
    return inputSpecifications.reduce((result, spec) => {
      let groupName = inputGroupName(spec.source);
      result[groupName] = result[groupName] || [];
      result[groupName].push(bindInputSpecification(spec,
        workflow,
        state.status,
        state.session,
        state.data,
        state.bufferedInputValues,
      dispatch));
      return result;
    }, {} as {[k: string]: BoundInputSpecification[]});
  }, [dispatch, inputSpecifications, state.bufferedInputValues, state.data, state.session, state.status, workflow]);

  return (
    <div className='primary-inputs'>
      {Object.keys(groupedInputs)
        .sort()
        .map((groupName, index) => {
          return (
            <InputGroup
              groupName={groupName}
              key={index}
              inputs={groupedInputs[groupName]}
            />
          );
        })}
    </div>
  );
}
