import React from 'react';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {WithInputParams} from './input_types';
import {maybeOfNullable, some, withDefault} from "../../../../selectors/maybe";
import {flattenStringOptions, StringOptions} from "./monoids";
import {useMemoized} from "../../../../selectors/workflow_selectors";
import { useSetsDefault } from './useSetsDefault';

export default function SelectAutocompleteInput({data, onChange, ...props}: WithInputParams<{}, string | null, StringOptions>) {
  const options = useMemoized(flattenStringOptions, data);
  const value = useSetsDefault("", props.value, onChange);

  return (
    <DropdownAutocomplete
      onSelect={(e: string | null) => {
        onChange(maybeOfNullable(e));
      }}
      list={options}
      value={value}
    />
  );
};
