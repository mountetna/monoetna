import React from 'react';
import {WithInputParams} from '../input_types';
import {some} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../useSetsDefault';
import {selectDefaultString} from '../monoids';
import { TextField } from '@material-ui/core';

export default function StringInput({onChange, label, minWidth, data, ...props}: WithInputParams<{label?: string, minWidth?: number}, string, string>) {
  const value = useSetsDefault(selectDefaultString(data), props.value.value, onChange);
  
  return (
    <div style={{paddingTop: label ? 8 : 0}}>
      <TextField
        value={value}
        multiline
        label={label}
        InputLabelProps={{ shrink: true }}
        onChange={(event) => onChange({value: some(event.target.value)})}
        size="small"
        style={{minWidth: minWidth || 200}}
      />
    </div>
  );
}