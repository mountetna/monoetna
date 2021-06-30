import React, {useCallback, useContext, useEffect, useMemo} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import {allWorkflowPrimaryInputSources, inputGroupName, useMemoized} from '../../../selectors/workflow_selectors';
import InputGroup from './input_group';
import {InputSpecification} from "../user_interactions/inputs/input_types";
import {useWorkflow} from "../../../contexts/workflow_context";
import {useCancelledOnDismount} from "etna-js/utils/cancellable";
import {Maybe, some, withDefault} from "../../../selectors/maybe";

export default function PrimaryInputs() {
  const {state, dispatch, startInputChange, onInputChange} = useContext(VulcanContext);
  const {session, bufferedInputValues} = state;
  const {workflow} = useWorkflow();

  const context = useCancelledOnDismount();

  const onStartInputChange = useCallback(() => {
    startInputChange(null, context);
  }, [context, startInputChange])

  const handlePrimaryInputChange = useCallback((inputName: string, val: any) => {
    onInputChange(null, {...bufferedInputValues, [inputName]: some(val)});
  }, [bufferedInputValues, onInputChange]);

  // Ensure defaults are set.
  useEffect(() => {
    let updatedInputs = bufferedInputValues;
    Object.keys(workflow.inputs).forEach(inputName => {
      if (!(inputName in bufferedInputValues)) {
        if (updatedInputs === bufferedInputValues) updatedInputs = {...updatedInputs};
        updatedInputs[inputName] = some(workflow.inputs[inputName].default);
      }
    })

    if (updatedInputs !== bufferedInputValues) {
      onInputChange(null, updatedInputs);
    }
  }, [bufferedInputValues, onInputChange, workflow.inputs])

  let groupedInputs = Object.entries(workflow.inputs).reduce((result, [inputName, input]) => {
    const inputSpecification: InputSpecification = {
      ...input,
      label: input.label || inputName,
      value: withDefault(bufferedInputValues[inputName], input.default),
      onChange: (v: Maybe<any>) => handlePrimaryInputChange(inputName, v),
    };

    let groupName = inputGroupName(inputName);
    result[groupName] = result[groupName] || [];
    result[groupName].push(inputSpecification);
    return result;
  }, {} as {[k: string]: InputSpecification[]});

  return (
    <div className='primary-inputs' onFocus={}>
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
