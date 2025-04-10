// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import {WithInputParams} from '../../input_types';
import {maybeOfNullable, some} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import NestedDropdownPiece from '../pieces/nested_dropdown_piece';
import { pullRecommendationIntoLabel } from './dropdown';
import { useMemo } from 'react';
import { nestedOptionSet } from '../pieces/utils';

export default function NestedDropdownInput({ label, data, onChange, ...props }: WithInputParams<{label?: string}, string|null, nestedOptionSet>) {
  const picked: string | null = useSetsDefault(null, props.value, onChange, 'picked');
  const allOptions = data.nestedOptions;
  const disp_label = useMemo(() => {
    return pullRecommendationIntoLabel(data, label);
  }, [data, label]);

  return NestedDropdownPiece(
    `nested-dropdown-${label}`,
    (v, k) => {onChange({picked: maybeOfNullable(v)})},
    picked,
    disp_label,
    allOptions,
    true
  )
};
