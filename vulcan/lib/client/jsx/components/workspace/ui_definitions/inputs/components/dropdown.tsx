import React, {useMemo} from 'react';
import {DataEnvelope, WithInputParams} from '../../input_types';
import {maybeOfNullable, some, withDefault} from '../../../../../selectors/maybe';
import {StringOptions} from '../../monoids';
import {useSetsDefault} from '../../useSetsDefault';
import { DropdownPieceRct } from '../pieces/dropdown_piece';

export function pullRecommendationIntoLabel<T extends DataEnvelope<any>>(
  data: T | null | undefined, label: string | undefined
): string | undefined {
  // Make the recommendation part of the autocomplete into the 'label' of renderInput
  if (data != null) {
    let data_use = {...data};
    let suggestion: string | undefined = undefined;
    if (Object.keys(data).includes('recommendation')) {
      // Allow multi-recommendation?
      const rec = Array.isArray(data['recommendation'])
        ? data['recommendation'].join(', ')
        : data['recommendation'];
      suggestion = rec == null || rec == 'null' ? undefined : 'Recommendation: ' + rec;
      delete data_use['recommendation'];
    }
    if (!!label) {
      return !!suggestion ? `${label}: ${suggestion}` : label
    }
    return suggestion;
  }
  return label;
}

export default function DropdownInput({
  data,
  defaultValue,
  label = '',
  minWidth,
  maxOptions = 200,
  disableClearable = true,
  disabled = false,
  showError,
  onChange,
  ...props
}: WithInputParams<
  {
    label?: string;
    minWidth?: number;
    maxOptions?: number;
    disableClearable?: boolean;
    disabled?: boolean;
  },
  string | null,
  StringOptions
>) {
  /*
  Creates a searchable dropdown selection input box from concatenated values of the 'data' hash.
  Special Case: If any data key is "recommendation", a line of text will display the values of this recommendation to the user.
  */
  const options_in: string[] = data.options;
  if (!options_in) showError('No data options given')

  const value = useSetsDefault(defaultValue || null, props.value, onChange, 'picked');
  const disp_label = useMemo(() => {
    return pullRecommendationIntoLabel(data, label);
  }, [data, label]);

  return <DropdownPieceRct
    name={`dropdown-${label}`}
    changeFxn={(v, k) => {onChange({picked: maybeOfNullable(v)})}}
    value={value}
    label={disp_label}
    options_in={options_in}
    minWidth={minWidth}
    disabled={disabled}
  />
}

