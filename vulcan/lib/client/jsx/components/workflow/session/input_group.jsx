import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import UserInput from '../user_interactions/inputs/user_input';
import {inputGroupName} from '../../../utils/workflow';
import {
  completedUiStepsSelector,
  nextUiStepsSelector
} from '../../../selectors/workflow';

export default function InputGroup({inputs, onChange}) {
  const context = useContext(VulcanContext);
  let {session} = context;

  const [open, setOpen] = useState(true);
  let groupName = inputGroupName(inputs[0]);

  function toggleInputs() {
    setOpen(!open);
  }

  useEffect(() => {
    if (session && session.inputs) {
      // Automatically collapse the input group if there are
      //   pending or completed UI input steps.
      let uiSteps = completedUiStepsSelector(context);
      let nextUiSteps = nextUiStepsSelector(context);

      if (uiSteps.length > 0 || nextUiSteps.length > 0) setOpen(false);
    }
  }, [session]);

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
