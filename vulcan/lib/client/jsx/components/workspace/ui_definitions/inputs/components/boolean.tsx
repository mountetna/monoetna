import React from 'react';
import {WithInputParams} from '../../input_types';
import {some} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { FormControlLabelProps } from '@material-ui/core';
import { CheckboxPieceRct } from '../pieces/user_input_pieces';

export default function BooleanInput({onChange, label, defaultValue, labelPlacement='end', disabled=false, data, ...props}: WithInputParams<{label?: string, labelPlacement?: FormControlLabelProps['labelPlacement'], disabled?: boolean}, boolean, boolean>) {
  const value = useSetsDefault(defaultValue || false, props.value, onChange, 'value');
  return <CheckboxPieceRct
    name={`boolean-${label}` || 'boolean-input'}
    changeFxn={(v, k?) => onChange({value: some(v)})}
    value={value}
    label={label || data.label || 'unlabeled-checkbox'}
    disabled={disabled}
    labelPlacement={labelPlacement}
  />
}