import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import ListInput from 'etna-js/components/inputs/list_input';
import SelectInput from 'etna-js/components/inputs/select_input';
import Dropdown from 'etna-js/components/inputs/dropdown';
import {
  IntegerInput,
  FloatInput
} from 'etna-js/components/inputs/numeric_input';
import SlowTextInput from 'etna-js/components/inputs/slow_text_input';
import Toggle from 'etna-js/components/inputs/toggle';

import {TYPE} from '../../../models/steps';

export default function PrimaryInputs() {
  const {workflow, session, setInputs} = useContext(VulcanContext);

  if (!workflow || !workflow.inputs || !session) return null;

  function handleInputChange(inputName, value) {
    let userInputs = {...session.inputs};
    userInputs[inputName] = value;
    setInputs(userInputs);
  }

  function wrapInput(inputComponent, input, inputName) {
    return (
      <div className='view_item'>
        <div className='item_name'>{input.label || inputName}</div>
        <div className='item_view'>{inputComponent}</div>
      </div>
    );
  }

  let components = Object.keys(workflow.inputs).map((inputName, index) => {
    let input = workflow.inputs[inputName];
    switch (input.type) {
      case TYPE.INTEGER:
        return wrapInput(
          <IntegerInput
            key={index}
            defaultValue={input.default}
            onChange={(e) => {
              handleInputChange(inputName, e);
            }}
          ></IntegerInput>,
          input,
          inputName
        );
      case TYPE.FLOAT:
        return wrapInput(
          <FloatInput
            key={index}
            defaultValue={input.default}
            onChange={(e) => {
              handleInputChange(inputName, e);
            }}
          ></FloatInput>,
          input,
          inputName
        );
        break;
      case TYPE.BOOL:
        return wrapInput(
          <input
            key={index}
            type='checkbox'
            className='text_box'
            onChange={(e) => {
              handleInputChange(inputName, e);
            }}
            defaultChecked={input.default}
          />,
          input,
          inputName
        );
        break;
      default:
        return wrapInput(
          <SlowTextInput key={index}></SlowTextInput>,
          input,
          inputName
        );
        break;
    }
  });

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
