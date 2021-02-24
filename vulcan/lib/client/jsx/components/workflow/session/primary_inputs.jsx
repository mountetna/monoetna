import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import {wrapEditableInputs} from '../../../utils/workflow';

export default function PrimaryInputs() {
  const {workflow, session, setInputs} = useContext(VulcanContext);

  if (!workflow || !workflow.inputs || !session) return null;

  function handleInputChange(inputName, value) {
    let userInputs = {...session.inputs};
    userInputs[inputName] = value;
    setInputs(userInputs);
  }

  let components = wrapEditableInputs(workflow.inputs, handleInputChange);

  return (
    <div className='primary-inputs inputs-pane'>
      <div class='title'>Inputs</div>
      <div className='primary-inputs-container items'>
        {components.map((comp) => {
          return comp;
        })}
      </div>
    </div>
  );
}
