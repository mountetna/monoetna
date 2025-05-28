import React from 'react';
import {WithInputParams} from '../../input_types';
import {maybeOfNullable, withDefault} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { NumberPieceRct } from '../pieces/number_pieces';

export function NumberInput(
  {
    onChange,
    label,
    minWidth,
    data,
    defaultValue,
    ...props
  }: WithInputParams<
    {minWidth?: number},
    number | null,
    number | null
  >,
  block_decimal = false
) {
  const value = useSetsDefault(defaultValue as number || 0, props.value, onChange, 'value'); // Had to replace selectDefaultNumber(data) with 0 due how the component can be cleared.  NaN was not an option because of cross-language conversion.
  return <NumberPieceRct
    name='number-input'
    changeFxn={(val: number | null) => {onChange({value: maybeOfNullable(val)})}}
    value={value}
    label={label}
    minWidth={minWidth}
    fullWidth={true}
    integer={block_decimal}
  />
}

export default function FloatInput({
  onChange,
  label,
  minWidth,
  data,
  defaultValue,
  ...props
}: WithInputParams<
  {minWidth?: number},
  number | null,
  number | null
>) {
  return NumberInput({onChange, label, minWidth, data, defaultValue, ...props});
}
