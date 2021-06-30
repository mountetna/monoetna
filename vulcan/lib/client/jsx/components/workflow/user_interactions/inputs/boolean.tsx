import React from 'react';
import {InputBackendComponent, WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";

export default function BooleanInput({input, onChange}: WithInputParams<{}>) {
  return (
    <input
      type='checkbox'
      className='text_box'
      onChange={(e) => {
        onChange(some(e.target.checked));
      }}
      checked={!!withDefault(input.value, false)}
    />
  );
}