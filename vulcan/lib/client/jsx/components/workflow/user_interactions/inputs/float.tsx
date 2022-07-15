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
  const [inputText, setInputText] = useState(from_num(value))
  const [sentValue, setSentValue] = useState(some(value))
  const [hasError, setHasError] = useState(sentValue==null || sentValue[0]==null)
  const onNewEntry = useCallback(
    (text: string, onlyInputState: boolean = false) => {
      const parser = ( block_decimal ) ? parseIntBetter : Number
      const parsed = parser(text)
      const newVal = (isNaN(parsed)) ? some(null) : some(parsed)
      setInputText(text)
      setSentValue(newVal)
      if (!onlyInputState) onChange(newVal)
    },
    [onChange]) // todo: better typing for event

  // Catch & show non-user value updates.
  useEffect(() => {
    if (props.value != sentValue) {
      onNewEntry(from_num(withDefault(props.value,null)), true)
    }
    setHasError(sentValue==null || sentValue[0]==null)
  }, [props.value])
  
  return (
    <div style={{paddingTop: label ? 8 : 0}}>
      <TextField
        value={inputText}
        label={label}
        error={hasError}
        onChange={(event) => onNewEntry(event.target.value)}
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