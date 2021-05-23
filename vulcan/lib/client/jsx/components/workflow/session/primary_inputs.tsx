import React, {useCallback, useContext, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {inputGroupName} from '../../../selectors/workflow_selectors';
import InputGroup from './input_group';
import {patchInputs} from "../../../actions/vulcan";
import {InputSpecification} from "../user_interactions/inputs/input_types";

export default function PrimaryInputs() {
  const {state, dispatch} = useContext(VulcanContext);
  const {workflow, session} = state;

  if (!workflow) return null;
  if (!workflow.inputs) return null;

  const handleInputChange = useCallback((inputName: string, val: any) => {
    dispatch(patchInputs({
      [inputName]: val,
    }));
  }, []);

  let primaryInputs: InputSpecification[] = useMemo(() => {
    return Object.keys(workflow.inputs).map(name => ({
        ...workflow.inputs[name],
          name,
          label: workflow.inputs[name].label || name,
          default: [
              session.inputs[name],
              workflow.inputs[name].default
          ].find(a => a != null)
    }))
  }, [session]);

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
