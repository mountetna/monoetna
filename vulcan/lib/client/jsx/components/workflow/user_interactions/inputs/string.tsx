import React from 'react';
import {WithInputParams} from './input_types';
import {some} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import {selectDefaultString} from "./monoids";
import { TextField } from '@material-ui/core';

export default function StringInput({onChange, label, data, ...props}: WithInputParams<{label?: string}, string, string>) {
  const value = useSetsDefault(selectDefaultString(data), props.value, onChange);
  
  if (label) {
    return <TextField
      value={value}
      label={label}
      variant="outlined"
      onChange={(event) => onChange(some(event.target.value))}
      size="small"
      style={{paddingTop:6}}
    />;
  }
  return <TextField
    value={value}
    variant="outlined"
    onChange={(event) => onChange(some(event.target.value))}
    size="small"
  />;
}