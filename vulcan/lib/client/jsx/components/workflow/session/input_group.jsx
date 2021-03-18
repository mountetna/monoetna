import React, {useState} from 'react';

import UserInput from '../user_interactions/inputs/user_input';
import {inputGroupName} from '../../../utils/workflow';

export default function InputGroup({inputs, onChange}) {
  const [open, setOpen] = useState(true);
  let groupName = inputGroupName(inputs[0]);

  function toggleInputs() {
    setOpen(!open);
  }

  return (
    <div className='inputs-pane'>
      <div
        onClick={toggleInputs}
        className={`inputs-pane-header ${open ? 'open' : 'closed'}`}
      >
        <div className='inputs-pane-toggle'></div>
        <div className='title'>{groupName}</div>
      </div>
      <div
        className={`primary-inputs-container items ${open ? 'open' : 'closed'}`}
      >
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
