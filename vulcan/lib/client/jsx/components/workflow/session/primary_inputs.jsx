import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import UserInput from '../user_interactions/user_input';

function InputGroup({inputs, onChange}) {
  let groupName = inputs[0].group || 'Inputs';

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
    if (!result.hasOwnProperty(input.group)) result[input.group] = [];

    result[input.group].push(input);

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
              onClick={handleInputChange}
            ></InputGroup>
          );
        })}
    </div>
  );
}
