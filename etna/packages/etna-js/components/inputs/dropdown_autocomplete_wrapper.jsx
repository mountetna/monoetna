// Wraps the DropdownAutocomplete component so it will work with the
// list_input to act as a multi-select.
import React, {useCallback} from 'react';

import DropdownAutocomplete from './dropdown_autocomplete';

export default function DropdownAutocompleteInput({
  onChange,
  list,
  onBlur,
  ...otherProps
}) {
  const onChangeWrapper = useCallback((value) => {
    onChange(value);
    onBlur();
  });

  return (
    <DropdownAutocomplete
      {...otherProps}
      list={list}
      onSelect={onChangeWrapper}
    />
  );
}
