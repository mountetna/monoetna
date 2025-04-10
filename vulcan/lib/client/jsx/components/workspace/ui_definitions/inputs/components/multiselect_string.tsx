import { MultiselectPiece } from '../pieces/user_input_pieces';
import {WithInputParams} from '../../input_types';
import {mapSome, maybeOfNullable, some, withDefault} from '../../../../../selectors/maybe';
import {StringOptions} from '../../monoids';

export default function MultiselectStringInput({onChange, data, label, defaultValue, ...props}: WithInputParams<
  {label?: string}, string[], string[]>) {
  const options = data.options;
  // const picked = useSetsDefault(defaultValue as string[] || [], props.value, onChange, 'picked'); // Had to replace selectDefaultNumber(data) with 0 due how the component can be cleared.  NaN was not an option because of cross-language conversion.
  const picked = withDefault(mapSome(props.value.picked, inner => Array.isArray(inner) ? inner : [inner]), []);

  return MultiselectPiece(
    label || 'multiselect',
    (val: string[] | null) => {onChange({picked: maybeOfNullable(val)})},
    picked,
    label || 'Selections',
    options,
    options,
    [],
    25
  );
}
