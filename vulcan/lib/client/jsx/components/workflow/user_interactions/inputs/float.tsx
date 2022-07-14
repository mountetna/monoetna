import React, {useCallback, useEffect, useState} from 'react';
import {WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";
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
  const value = useSetsDefault(0, props.value, onChange); // Had to replace selectDefaultNumber(data) with 0 due how the component can be cleared.  NaN was not an option because of cross-language conversion.
  const [inputState, setInputState] = useState({text: from_num(value), hasError: false})
  const onNewFloat = useCallback(
    (value: string, onlyInputState: boolean = false) => {
      const parser = ( block_decimal ) ? parseIntBetter : Number
      const parsed = parser(value)
      if (isNaN(parsed)) {
        setInputState({text: value, hasError: true})
        if (!onlyInputState) onChange([null])
      } else {
        setInputState({text: value, hasError: false})
        if (!onlyInputState) onChange(some(parsed))
      }
    },
    [onChange]) // todo: better typing for event

  // Catch & show non-user value updates.
  useEffect(() => {
    onNewFloat(from_num(withDefault(props.value,null)), true)
  }, [props.value])
  
  return (
    <div style={{paddingTop: label ? 8 : 0}}>
      <TextField
        value={inputState.text}
        label={label}
        error={inputState.hasError}
        onChange={(event) => onNewFloat(event.target.value)}
        size="small"
        style={{minWidth: minWidth || 200}}
      />
    </div>
  )
};

function parseIntBetter(s: string){
  const parsed = parseInt(s, 10);
  return parsed + "" == s ? parsed : NaN;
}

function from_num(num: number | null) {
  return num == null ? "" : num.toString()
}