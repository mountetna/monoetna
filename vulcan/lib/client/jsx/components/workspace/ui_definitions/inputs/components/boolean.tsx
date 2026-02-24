import React from 'react';
import {WithInputParams} from '../../input_types';
import {some} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { FormControlLabelProps } from '@material-ui/core';
import { CheckboxPieceRct } from '../pieces/checkbox_piece';

export default function BooleanInput({
  onChange,
  label,
  defaultValue,
  labelPlacement='end',
  disabled=false,
  data,
  ...props
}: WithInputParams<
  {
    label?: string,
    labelPlacement?: FormControlLabelProps['labelPlacement'],
    disabled?: boolean
  },
  boolean, 
  any
>) {
  const value = useSetsDefault(defaultValue || false, props.value, onChange, 'value');
  const label_use = !!label ? label :
    !! data && ('label' in data) ? data['label'] as string :
    'unlabeled-checkbox';
  return <CheckboxPieceRct
    name={`boolean-${label}` || 'boolean-input'}
    changeFxn={(v, k?) => onChange({value: some(v)})}
    value={value}
    label={label_use}
    disabled={disabled}
    labelPlacement={labelPlacement}
  />
}