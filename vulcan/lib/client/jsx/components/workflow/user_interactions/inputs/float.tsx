import React, {useCallback} from 'react';
import {WithInputParams} from './input_types';
import {some} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import {selectDefaultNumber} from "./monoids";
import { TextField } from '@material-ui/core';

export function FloatInput({onChange, data, ...props}: WithInputParams<{}, number, number>) {
  const onNewFloat = useCallback((event: any) => onChange(some(event.target.value)), [onChange]) // better typing for event
  const value = useSetsDefault(selectDefaultNumber(data), props.value, onChange);

  return <TextField
    value={value}
    type="number"
    variant="outlined"
    onChange={onNewFloat}
    size="small"
  />;
}

export default FloatInput;
