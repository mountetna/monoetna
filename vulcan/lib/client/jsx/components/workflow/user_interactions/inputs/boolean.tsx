import React, {ChangeEvent, useCallback, useEffect} from 'react';
import {WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";

export default function BooleanInput({value, onChange, label, d = false}: WithInputParams<{label?: string, d?: boolean}, boolean>) {
  const onCheck = useCallback((e: ChangeEvent<HTMLInputElement>) => onChange(some(e.target.checked)), [onChange]);

  useEffect(() => {
    if (!value) onChange(some(d))
  }, [onChange, value, d])

  const inner = <input
      type='checkbox'
      className='text_box'
      onChange={onCheck}
      checked={withDefault(value, false)}
    />

  if (label) {
    return <label className='checkbox-input-option'>
      {inner}
    </label>;
  }

  return inner;
}