import React from 'react';
import {WithInputParams} from '../../input_types';
import {mapSome, maybeOfNullable, some, withDefault} from '../../../../../selectors/maybe';
import {StringOptions} from '../../monoids';
import { MultiselectStringPieceRct } from '../pieces/multiselect_string_piece';

export default function MultiselectStringInput({onChange, data, label, defaultValue, ...props}: WithInputParams<
  {}, string[], string[]>) {
  if (!data || !('options' in data)) {
    props.showError('required input data missing')
    return null
  }
  const options: string[] = data.options;
  const picked: string[] = withDefault(mapSome(props.value.picked, (inner: string | string[]) => Array.isArray(inner) ? inner : [inner]), []);

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
