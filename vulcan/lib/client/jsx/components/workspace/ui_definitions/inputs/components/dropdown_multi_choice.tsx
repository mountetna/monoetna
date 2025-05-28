import React, {useCallback, useEffect, useMemo, useState} from 'react';
import {WithInputParams} from '../../input_types';
import {useSetsDefault} from '../../useSetsDefault';
import { pullRecommendationIntoLabel } from './dropdown'
import { DropdownMultiChoicePieceRct } from '../pieces/dropdown_multi_choice_piece';
import { maybeOfNullable } from '../../../../../../../../../etna/packages/etna-js/selectors/maybe';

export default function DropdownMultiChoiceInput({
  data,
  label,
  defaultValue,
  onChange,
  ...props
}: WithInputParams<
  {},
  string[],
  string[]
>) {
  if (!!defaultValue && typeof defaultValue === 'string') {
    defaultValue = [defaultValue];
  }
  const options_in: string[] = data.options;

  const value = useSetsDefault(defaultValue || null, props.value, onChange, 'picked');
  const disp_label = useMemo(() => {
    return pullRecommendationIntoLabel(data, label);
  }, [data, label]);

  return (
    <DropdownMultiChoicePieceRct
      name={label || `dropdown-multi-choice-input`}
      changeFxn={(v: string[] | null, k?: string) => onChange({picked: maybeOfNullable(v)})}
      value={value}
      label={disp_label}
      options_in={options_in}
      disableClearable={false}
    />
  );
}
