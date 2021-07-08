import React, {ChangeEvent, useCallback} from 'react';
import {InputBackendComponent, WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";

export default function BooleanInput({value, onChange, label}: WithInputParams<{label?: string}>) {
  const onCheck = useCallback((e: ChangeEvent<HTMLInputElement>) => onChange(some(e.target.checked)), [onChange]);

  const inner = <input
      type='checkbox'
      className='text_box'
      onChange={onCheck}
      checked={!!withDefault(value, false)}
    />

  if (label) {
    return <label className='checkbox-input-option'>
      {inner}
    </label>;
  }

  return inner;
}