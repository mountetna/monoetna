import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import {wrapEditableInputs} from '../../../utils/workflow';

function InputGroup({inputs, onClick}) {
  let groupName = inputs[0].group || 'Inputs';

  let components = wrapEditableInputs(inputs, onClick);

  return (
    <div className='inputs-pane'>
      <div className='title'>{groupName}</div>
      <div className='primary-inputs-container items'>
        {components.map((comp) => {
          return comp;
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
          name: inputName,
          type: workflowInput.type,
          label: workflowInput.label || inputName,
          default:
            (session.inputs && session.inputs[inputName]) ||
            workflowInput.default ||
            null,
          group: workflowInput.group
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
