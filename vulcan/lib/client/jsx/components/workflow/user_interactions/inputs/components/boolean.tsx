import React, {ChangeEvent, useCallback, useEffect} from 'react';
import {WithInputParams} from '../input_types';
import {some, withDefault} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../useSetsDefault';
import {selectDefaultBoolean} from '../monoids';
import { Checkbox, FormControlLabel, FormControlLabelProps } from '@material-ui/core';

export default function BooleanInput({onChange, label, labelPlacement='end', disabled=false, data, ...props}: WithInputParams<{label?: string, labelPlacement?: FormControlLabelProps['labelPlacement'], disabled?: boolean}, boolean, boolean>) {
  const value = useSetsDefault(selectDefaultBoolean(data), props.value.value, onChange);
  const onCheck = useCallback((e: ChangeEvent<HTMLInputElement>) => onChange({value: some(!value)}), [onChange, value]);

  const inner = <Checkbox
      onChange={onCheck}
      checked={value}
      inputProps={{ 'aria-label': 'controlled' }}
      disabled={disabled}
    />;

  if (label) {
    return (
      <FormControlLabel control={inner} label={label} labelPlacement={labelPlacement}/>
    );
  }

  return inner;
}