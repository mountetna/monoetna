import React from 'react';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {InputBackendComponent} from "./input_types";

const SelectAutocompleteInput: InputBackendComponent = ({input, onChange}) => {
  if (!input || !onChange) return null;

    const options: any[] = Object.values(input.data || {}).reduce((acc, n) => {
        if (Array.isArray(n)) return acc.concat(n);
        return acc.concat([n]);
    }, [input.data]);

  return (
    <DropdownAutocomplete
      onSelect={(e: any) => {
        onChange(input.name, e);
      }}
      list={options}
      defaultValue={input.default || null}
    />
  );
}

export default SelectAutocompleteInput;
