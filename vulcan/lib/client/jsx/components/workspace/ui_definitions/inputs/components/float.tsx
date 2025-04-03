import {WithInputParams} from '../../input_types';
import {maybeOfNullable, withDefault} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { FloatPiece, IntegerPiece } from '../pieces/number_pieces';

export function NumberInput(
  {
    onChange,
    label,
    minWidth,
    data,
    defaultValue,
    ...props
  }: WithInputParams<
    {label?: string; minWidth?: number},
    number | null,
    number | null
  >,
  block_decimal = false
) {
  const value = useSetsDefault(defaultValue as number || 0, props.value, onChange, 'value'); // Had to replace selectDefaultNumber(data) with 0 due how the component can be cleared.  NaN was not an option because of cross-language conversion.
  const UI = block_decimal ? IntegerPiece : FloatPiece;
  return UI(
    'input', (val: number | null) => {onChange({value: maybeOfNullable(val)})},
    value,
    label, minWidth
  )
}

export default function FloatInput({
  onChange,
  label,
  minWidth,
  data,
  defaultValue,
  ...props
}: WithInputParams<
  {label?: string; minWidth?: number},
  number | null,
  number | null
>) {
  return NumberInput({onChange, label, minWidth, data, defaultValue, ...props});
}
