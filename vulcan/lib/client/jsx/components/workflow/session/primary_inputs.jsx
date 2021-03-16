import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import {wrapEditableInputs} from '../../../utils/workflow';

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

        result[inputName] = {
          type: workflowInput.type,
          label: workflowInput.label || inputName,
          default: session.inputs[inputName] || workflowInput.default || null
        };
        return result;
      },
      {}
    );

    return mergedInputs;
  }

  let components = wrapEditableInputs(
    mergeSessionDefaultInputs(),
    handleInputChange
  );

  return (
    <div className='primary-inputs inputs-pane'>
      <div className='title'>Inputs</div>
      <div className='primary-inputs-container items'>
        {components.map((comp) => {
          return comp;
        })}
      </div>
    </div>
  );
}
