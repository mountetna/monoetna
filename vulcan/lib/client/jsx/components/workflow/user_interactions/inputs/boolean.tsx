import React, {ChangeEvent, useCallback, useEffect} from 'react';
import {WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import {selectDefaultBoolean} from "./monoids";
import { Checkbox, FormControlLabel } from '@material-ui/core';

export default function BooleanInput({onChange, label, data, ...props}: WithInputParams<{label?: string}, boolean, boolean>) {
  const value = useSetsDefault(selectDefaultBoolean(data), props.value, onChange);
  const onCheck = useCallback((e: ChangeEvent<HTMLInputElement>) => onChange(some(!value)), [onChange, value]);

  const inner = <Checkbox
      onChange={onCheck}
      checked={value}
      inputProps={{ 'aria-label': 'controlled' }}
    />

  if (label) {
    return (
      <FormControlLabel control={inner} label={label}/>
    );
  }

  return inner;
}