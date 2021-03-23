import React from 'react';
import {InputBackendComponent} from "./types";


const BooleanInput: InputBackendComponent = ({input, onChange}) => {
  if (!input || !onChange) return null;

  return (
    <input
      type='checkbox'
      className='text_box'
      onChange={(e) => {
        onChange(input.name, e);
      }}
      defaultChecked={input.default}
    />
  );
}

export default BooleanInput;
