import React from 'react';
import {WithInputParams} from '../../input_types';
import {some} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import {selectDefaultString} from '../../monoids';
import { TextField } from '@material-ui/core';
import { stringPiece } from '../pieces/user_input_pieces';

export default function StringInput({onChange, label, minWidth, defaultValue, data, ...props}: WithInputParams<{label?: string, minWidth?: number}, string, string>) {
  const value = useSetsDefault(defaultValue || selectDefaultString(data), props.value, onChange, 'value');
  return stringPiece(
    label || 'string-input',
    (v, k?) => onChange({value: some(v)}),
    value,
    label,
    minWidth
  )
}
