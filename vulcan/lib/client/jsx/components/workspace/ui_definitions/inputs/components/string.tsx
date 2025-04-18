import React from 'react';
import {WithInputParams} from '../../input_types';
import {some} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { StringPieceRct } from '../pieces/string_piece';

export default function StringInput({onChange, label, minWidth, defaultValue, data, ...props}: WithInputParams<{label?: string, minWidth?: number}, string, string>) {
  const value = useSetsDefault(defaultValue || '', props.value, onChange, 'value');
  return <StringPieceRct
    name={label || 'string-input'}
    changeFxn={(v, k?) => onChange({value: some(v)})}
    value={value}
    label={label}
    minWidth={minWidth}
  />
}
