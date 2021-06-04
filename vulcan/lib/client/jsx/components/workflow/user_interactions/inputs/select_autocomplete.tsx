import React, {useMemo} from 'react';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {InputBackendComponent} from './input_types';

const SelectAutocompleteInput: InputBackendComponent = ({input, onChange}) => {
  const options: any[] = useMemo(
    () =>
      Object.values(input.data || {}).reduce((acc, n) => {
        if (Array.isArray(n)) return acc.concat(n);
        return acc;
      }, []),
    [input.data]
  );

  if (!input || !onChange) return null;

  return (
    <DropdownAutocomplete
      onSelect={(e: any) => {
        onChange(input.name, e);
      }}
      list={options}
      value={input.value || null}
    />
  );
};

export default SelectAutocompleteInput;
