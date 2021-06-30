import React, {useCallback, useContext, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {allWorkflowPrimaryInputSources, inputGroupName, useMemoized} from '../../../selectors/workflow_selectors';
import InputGroup from './input_group';
import {InputSpecification} from "../user_interactions/inputs/input_types";
import {useWorkflow} from "../../../contexts/workflow_context";

export default function PrimaryInputs() {
  const {state, dispatch} = useContext(VulcanContext);
  const {session} = state;
  const {workflow} = useWorkflow();

  const handleInputChange = useCallback((inputName: string, val: any) => {
    dispatch(patchInputs({
      [inputName]: val,
    }));
  }, [dispatch]);

  let primaryInputs: InputSpecification[] = useMemo(() => {
    return Object.keys(workflow.inputs).map(name => ({
        ...workflow.inputs[name],
          label: workflow.inputs[name].label || name,
          value: [
              session.inputs[name],
              workflow.inputs[name].default
          ].find(a => a != null)
    }))
  }, [session.inputs, workflow.inputs]);

  let groupedInputs = primaryInputs.reduce((result, input) => {
    let groupName = inputGroupName(input.name);
    result[groupName] = result[groupName] || [];
    result[groupName].push(input);
    return result;
  }, {} as {[k: string]: InputSpecification[]});

  return (
    <div className='primary-inputs'>
      {Object.keys(groupedInputs)
        .sort()
        .map((groupName, index) => {
          return (
            <InputGroup
              key={index}
              inputs={groupedInputs[groupName]}
              onChange={handleInputChange}
            />
          );
        })}
    </div>
  );
}
