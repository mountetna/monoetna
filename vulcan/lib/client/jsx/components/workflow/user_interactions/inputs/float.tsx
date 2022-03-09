import React, {useCallback} from 'react';
import {WithInputParams} from './input_types';
import {some} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import { TextField } from '@material-ui/core';

export default function FloatInput({
  onChange, label, minWidth, data, ...props
}: WithInputParams<{label?: string, minWidth?: number}, number | null, number | null>) {
    return NumberInput({onChange, label, minWidth, data, ...props})
}

export function NumberInput(
  {onChange, label, minWidth, data, ...props}: WithInputParams<{label?: string, minWidth?: number}, number | null, number | null>,
  block_decimal=false
  ) {
  const onNewFloat = useCallback(
    (event: any) => {
      let next_value = event.target.value
      if ( block_decimal && /\./.test(event.target.value) ) {
        const before_decimal = (/(.*)\./.exec(event.target.value) as string[])[1]
        next_value = /^[+-]?$/.test(before_decimal) ? "0" : before_decimal
      }
      onChange(some(to_num(next_value)))
    },
    [onChange]) // todo: better typing for event
  const value = useSetsDefault(0, props.value, onChange); // Had to replace selectDefaultNumber(data) with 0 due how the component can be cleared.  NaN was not an option because of cross-language conversion.
  
  return <TextField
    value={from_num(value)}
    label={label}
    type="number"
    onChange={onNewFloat}
    size="small"
    style={{minWidth: minWidth || 200}}
  />;
};

// Being a string input at heart, the MaterialUI TextField can be emptied with a backspace.
function to_num(str:string) {
  return str == "" ? null : Number(str)
}

function from_num(num: number | null) {
  return num == null ? "" : num.toString()
}