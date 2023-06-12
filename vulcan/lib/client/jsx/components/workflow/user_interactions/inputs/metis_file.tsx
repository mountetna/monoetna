import React, {ChangeEvent, useCallback, useEffect} from 'react';
import {WithInputParams} from './input_types';
import {some, withDefault} from '../../../../selectors/maybe';
import { PickFileOrFolder, PickBucket } from 'etna-js/components/metis_exploration';
import {useSetsDefault} from './useSetsDefault';
import {selectDefaultBoolean} from './monoids';
import { Checkbox, FormControlLabel, FormControlLabelProps } from '@material-ui/core';

export default function FileInput({onChange, label, data, ...props}: WithInputParams<{label?: string}, boolean, boolean>) {
  const value = useSetsDefault('', props.value, onChange);
  const onCheck = useCallback((e: ChangeEvent<HTMLInputElement>) => onChange(some(!value)), [onChange, value]);

  const inner = <Checkbox
      onChange={onCheck}
      checked={value}
      inputProps={{ 'aria-label': 'controlled' }}
    />;

  if (label) {
    return (
      <FormControlLabel control={inner} label={label} labelPlacement={labelPlacement}/>
    );
  }

  return inner;
}