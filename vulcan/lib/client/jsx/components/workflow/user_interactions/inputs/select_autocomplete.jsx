import React from 'react';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';

export default function SelectAutocompleteInput({input, onChange}) {
  if (!input || !onChange) return null;

  return (
    <DropdownAutocomplete
      onSelect={(e) => {
        onChange(input.name, e);
      }}
      list={input.options || []}
      defaultValue={input.default || null}
    ></DropdownAutocomplete>
  );
}
