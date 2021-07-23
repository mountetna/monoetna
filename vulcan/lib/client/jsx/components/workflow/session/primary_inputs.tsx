import React, {useCallback, useContext, useEffect, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {inputGroupName} from '../../../selectors/workflow_selectors';
import InputGroup from './input_group';
import {
  bindInputSpecification, BoundInputSpecification, DataEnvelope, getInputSpecifications
} from "../user_interactions/inputs/input_types";
import {useWorkflow} from "../../../contexts/workflow_context";
import {Maybe, maybeOfNullable} from "../../../selectors/maybe";
import {BufferedInputsContext, WithBufferedInputs} from "../../../contexts/input_state_management";

export default function PrimaryInputs() {
  const {commitSessionInputChanges, dispatch} = useContext(VulcanContext);

  return <WithBufferedInputs commitSessionInputChanges={commitSessionInputChanges} dispatch={dispatch} stepName={null}>
    <PrimaryInputsInner/>
  </WithBufferedInputs>
}


function PrimaryInputsInner() {
  const {state, dispatch} = useContext(VulcanContext);
  const {inputs, setInputs} = useContext(BufferedInputsContext);
  const {session} = state;
  const {workflow} = useWorkflow();

  // Ensure defaults are set.
  useEffect(() => {
    let withDefaults: DataEnvelope<Maybe<any>> = {};
    Object.keys(workflow.inputs).forEach(inputName => {
      if (!(inputName in session.inputs) && !(inputName in inputs)) {
        withDefaults[inputName] = maybeOfNullable(workflow.inputs[inputName].default);
      }
    })

    if (Object.keys(withDefaults).length > 0) {
      setInputs(inputs => ({...inputs, ...withDefaults}));
    }
  }, [inputs, session.inputs, setInputs, workflow.inputs])

  const inputSpecifications = useMemo(() =>
    getInputSpecifications(Object.entries(workflow.inputs), workflow), [workflow]);

  let groupedInputs = useMemo(() => {
    return inputSpecifications.reduce((result, spec) => {
      let groupName = inputGroupName(spec.source) || "Inputs";
      result[groupName] = result[groupName] || [];
      result[groupName].push(bindInputSpecification(spec,
        workflow,
        state.status,
        state.session,
        state.data,
        inputs,
        setInputs,
      ));
      return result;
    }, {} as {[k: string]: BoundInputSpecification[]});
  }, [inputSpecifications, inputs, setInputs, state.data, state.session, state.status, workflow]);

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