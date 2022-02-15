import React, {ChangeEvent, useCallback, useEffect} from 'react';
import {WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import {selectDefaultBoolean} from "./monoids";

export default function BooleanInput({onChange, label, data, ...props}: WithInputParams<{label?: string}, boolean, boolean>) {
  const value = useSetsDefault(selectDefaultBoolean(data), props.value, onChange);
  const onCheck = useCallback((e: ChangeEvent<HTMLInputElement>) => onChange(some(!value)), [onChange, value]);

  const inner = <input
      type='checkbox'
      className='text_box'
      onChange={onCheck}
      checked={value}
    />

  if (label) {
    return <label className='checkbox-input-option' style={{display: 'inline-flex'}}>
      {inner}
      {label}
    </label>;
  }

  return inner;
}