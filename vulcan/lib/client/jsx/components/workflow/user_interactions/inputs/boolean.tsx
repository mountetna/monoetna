import React from 'react';
import {InputBackendComponent} from './input_types';

const BooleanInput: InputBackendComponent = ({input, onChange}) => {
  if (!input || !onChange) return null;

  return (
    <input
      type='checkbox'
      className='text_box'
      onChange={(e) => {
        onChange(input.name, e.target.checked);
      }}
      checked={input.value || false}
    />
  );
};

export default BooleanInput;
