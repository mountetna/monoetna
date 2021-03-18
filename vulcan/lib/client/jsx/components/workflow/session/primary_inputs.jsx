import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import UserInput from '../user_interactions/inputs/user_input';
import {inputGroupName} from '../../../utils/workflow';

function InputGroup({inputs, onChange}) {
  let groupName = inputGroupName(inputs[0]);

  return (
    <div className='inputs-pane'>
      <div className='title'>{groupName}</div>
      <div className='primary-inputs-container items'>
        {inputs.map((input, index) => {
          return (
            <UserInput
              input={input}
              onChange={onChange}
              key={index}
            ></UserInput>
          );
        })}
      </div>
    </div>
  );
}

export default function PrimaryInputs() {
  const {workflow, session, setInputs} = useContext(VulcanContext);

  if (!workflow || !workflow.inputs) return null;

  function handleInputChange(inputName, value) {
    let userInputs = {...session.inputs};
    userInputs[inputName] = value;
    setInputs(userInputs);
  }

  function mergeSessionDefaultInputs() {
    // Merge the session and workflow default inputs,
    //   with any session inputs taking precedence.
    let mergedInputs = Object.keys(workflow.inputs).reduce(
      (result, inputName) => {
        let workflowInput = workflow.inputs[inputName];

        result.push({
          ...workflowInput,
          name: inputName,
          label: workflowInput.label || inputName,
          default:
            (session.inputs && session.inputs[inputName]) ||
            workflowInput.default ||
            null
        });

        return result;
      },
      []
    );

    return mergedInputs;
  }

  let primaryInputs = mergeSessionDefaultInputs();

  let groupedInputs = primaryInputs.reduce((result, input) => {
    let groupName = inputGroupName(input);

    if (!result.hasOwnProperty(groupName)) result[groupName] = [];

    result[groupName].push(input);

    return result;
  }, {});

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
            ></InputGroup>
          );
        })}
    </div>
  );
}
