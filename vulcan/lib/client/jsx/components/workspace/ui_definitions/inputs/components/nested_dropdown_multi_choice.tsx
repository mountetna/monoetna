// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useMemo} from 'react';
import {WithInputParams} from '../../input_types';
import {maybeOfNullable} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { nestedOptionSet } from '../pieces/utils';
import NestedDropdownMultiChoicePieceRct from '../pieces/nested_dropdown_multi_choice_piece';
import { pullRecommendationIntoLabel } from './dropdown';
import { NestedDropdownMultiChoiceAdvancedPieceRct } from '../pieces/nested_dropdown_multi_choice_advanced_piece';

function _NestedDropdownMultiChoiceInput({ letReorder, letBulkAdd, label, testIdAppend='', data, onChange, ...props }: WithInputParams<{
  letReorder: boolean;
  letBulkAdd: boolean;
  label?: string;
  testIdAppend?: string;
}, string[], nestedOptionSet>) {
  const value: string[] = useSetsDefault([] as string[], props.value, onChange, 'picked');
  const options_in = data.nestedOptions;
  const disp_label = useMemo(() => {
    return pullRecommendationIntoLabel(data, label);
  }, [data, label]);

  if (letReorder || letBulkAdd) {
    return <NestedDropdownMultiChoiceAdvancedPieceRct
      name={label || 'nested-dropdown-multi-choice-advanced-input'}
      changeFxn={(v,k) => {onChange({picked: maybeOfNullable(v)})}}
      value={value}
      label={disp_label}
      options_in={options_in}
      letBulkAdd={letBulkAdd}
      letReorder={letReorder}
    />
  }
  return <NestedDropdownMultiChoicePieceRct
    name={label || 'nested-dropdown-multi-choice-input'}
    changeFxn={(v,k) => {onChange({picked: maybeOfNullable(v)})}}
    value={value}
    label={disp_label}
    options_in={options_in}
  />
};

export function NestedDropdownMultiChoiceInput({ label, testIdAppend='', data, onChange, ...props }: WithInputParams<{
  label?: string
  testIdAppend?: string
}, string[], nestedOptionSet>) {
  return  _NestedDropdownMultiChoiceInput({ letReorder: false, letBulkAdd: false, label, testIdAppend, data, onChange, ...props });
};

export function NestedDropdownMultiChoiceBulkAddInput({ label, testIdAppend='', data, onChange, ...props }: WithInputParams<{
  label?: string
  testIdAppend?: string
}, string[], nestedOptionSet>) {
  return  _NestedDropdownMultiChoiceInput({ letReorder: false, letBulkAdd: true, label, testIdAppend, data, onChange, ...props });
};

export function NestedDropdownMultiChoiceReorderInput({ label, testIdAppend='', data, onChange, ...props }: WithInputParams<{
  label?: string
  testIdAppend?: string
}, string[], nestedOptionSet>) {
  return  _NestedDropdownMultiChoiceInput({ letReorder: true, letBulkAdd: false, label, testIdAppend, data, onChange, ...props });
};

export function NestedDropdownMultiChoiceAdvancedInput({ label, testIdAppend='', data, onChange, ...props }: WithInputParams<{
  label?: string
  testIdAppend?: string
}, string[], nestedOptionSet>) {
  return  _NestedDropdownMultiChoiceInput({ letReorder: true, letBulkAdd: true, label, testIdAppend, data, onChange, ...props });
};
