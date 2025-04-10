import React from 'react';
import {WithInputParams} from '../../input_types';
import {mapSome, maybeOfNullable, some, withDefault} from '../../../../../selectors/maybe';
import {StringOptions} from '../../monoids';
import { MultiselectStringPieceRct } from '../pieces/multiselect_string_piece';

export default function MultiselectStringInput({onChange, data, label, defaultValue, ...props}: WithInputParams<
  {label?: string}, string[], string[]>) {
  const options: string[] = data.options;
  // const picked = useSetsDefault(defaultValue as string[] || [], props.value, onChange, 'picked'); // Had to replace selectDefaultNumber(data) with 0 due how the component can be cleared.  NaN was not an option because of cross-language conversion.
  const picked = withDefault(mapSome(props.value.picked, inner => Array.isArray(inner) ? inner : [inner]), []);

  return <MultiselectStringPieceRct
    name={label || 'multiselect-input'}
    changeFxn={(val: string[] | null) => {onChange({picked: maybeOfNullable(val)})}}
    value={picked}
    label={label || 'Selections'}
    options={options}
    onAll={options}
    onClear={[]}
    maxItems={25}
  />;
}
