import React, {useMemo} from 'react';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";

export default function SelectAutocompleteInput({data, onChange, ...props}: WithInputParams<{}>) {
  const value = withDefault(props.value, null);
  const options: any[] = useMemo(
    () =>
      Object.values(data || {}).reduce((acc, n) => {
        if (Array.isArray(n)) return acc.concat(n);
        return acc;
      }, []),
    [data]
  );

  return (
    <DropdownAutocomplete
      onSelect={(e: string | null) => {
        onChange(some(e));
      }}
      list={options}
      value={value}
    />
  );
};
