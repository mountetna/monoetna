import React, {useState, useContext, useEffect, useCallback} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';

import UserInput from '../user_interactions/inputs/user_input';
import {
  inputGroupName,
  nextUiStepsSelector
} from '../../../selectors/workflow_selectors';
import {InputSpecification} from "../user_interactions/inputs/types";

export default function InputGroup({inputs, onChange}: {inputs: InputSpecification[], onChange: (name: string, val: any) => void}) {
  const {state} = useContext(VulcanContext);
  let {session} = state;
  const [open, setOpen] = useState(true);

  if (inputs.length === 0) return null;
  let groupName = inputGroupName(inputs[0].name);

  const toggleInputs = useCallback(() => setOpen(!open), [setOpen, open]);

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
      <div className='header-wrapper'>
        <div onClick={toggleInputs} className='inputs-pane-header'>
          <div className='title'>{groupName}</div>
        </div>
        <div className='filler'/>
      </div>
      <div
        className={`primary-inputs-container items sliding-panel vertical ${
          open ? 'open' : 'closed'
        }`}
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
