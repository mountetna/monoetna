// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React from 'react';
import {nestedOptionSet, WithInputParams} from '../../input_types';
import {maybeOfNullable, some} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { NestedDropdownPieceRct } from '../pieces/nested_dropdown_piece';
import { pullRecommendationIntoLabel } from './dropdown';
import { useMemo } from 'react';

export default function NestedDropdownInput({ label, data, onChange, ...props }: WithInputParams<{}, string|null, nestedOptionSet>) {
  const picked: string | null = useSetsDefault(null, props.value, onChange, 'picked');
  const allOptions = data.nestedOptions;
  const disp_label = useMemo(() => {
    return pullRecommendationIntoLabel(data, label);
  }, [data, label]);

  return <NestedDropdownPieceRct
    name={`nested-dropdown-${label}`}
    changeFxn={(v, k) => {onChange({picked: maybeOfNullable(v)})}}
    value={picked}
    label={disp_label}
    nestedOptions={allOptions}
  />
};
