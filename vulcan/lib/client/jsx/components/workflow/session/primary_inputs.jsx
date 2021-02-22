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

  let components = Object.keys(workflow.inputs).map((inputName, index) => {
    let input = workflow.inputs[inputName];
    switch (input.type) {
      case TYPE.INTEGER:
        return (
          <div className='view_item'>
            <div className='item_name'>{input.label || inputName}</div>
            <div className='item_view'>
              <IntegerInput
                key={index}
                defaultValue={input.default}
                onChange={(e) => {
                  handleInputChange(inputName, e);
                }}
              ></IntegerInput>
            </div>
          </div>
        );
        break;
      case TYPE.FLOAT:
        return <FloatInput key={index}></FloatInput>;
        break;
      case TYPE.BOOL:
        return <Toggle key={index}></Toggle>;
        break;
      default:
        return <SlowTextInput key={index}></SlowTextInput>;
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
