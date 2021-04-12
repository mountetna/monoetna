// Wraps the Dropdown component so it will work with the
// list_input to act as a multi-select.
import React from 'react';

import DropdownAutocomplete from './dropdown_autocomplete';

export default function DropdownAutocompleteInput({
  onChange,
  list,
  onBlur,
  ...otherProps
}) {
  function onChangeWrapper(value) {
    onChange(value);
    onBlur();
  }

  return (
    <DropdownAutocomplete
      {...otherProps}
      list={list}
      onSelect={onChangeWrapper}
    />
  );
}
