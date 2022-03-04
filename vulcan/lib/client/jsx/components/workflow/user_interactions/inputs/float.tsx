import React, {useCallback} from 'react';
import {WithInputParams} from './input_types';
import {some} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import {selectDefaultNumber} from "./monoids";
import { TextField } from '@material-ui/core';

export function FloatInput({onChange, label, minWidth, data, ...props}: WithInputParams<{label?: string, minWidth?: number}, number | null, number | null>) {
  const onNewFloat = useCallback((event: any) => onChange(some(to_num(event.target.value))), [onChange]) // todo: better typing for event
  const value = useSetsDefault(0, props.value, onChange); // Had to replace selectDefaultNumber(data) with 0 due how the component can be cleared.  NaN was not an option because of cross-language conversion.

  return <TextField
    value={from_num(value)}
    label={label}
    type="number"
    variant="outlined"
    onChange={onNewFloat}
    size="small"
    style={{minWidth: minWidth || 200}}
  />;
}

// Being a string input at heart, the MaterialUI TextField can be emptied with a backspace.
function to_num(str:string) {
  return str == "" ? null : Number(str)
}

function from_num(num: number | null) {
  return num == null ? "" : num.toString()
}

export default FloatInput;
