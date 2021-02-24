// Wraps the Dropdown component so it will work with the
// list_input to act as a multi-select.
import React from 'react';

import Dropdown from './dropdown';

export default function DropdownInput({onChange, list, onBlur, ...otherProps}) {
  function onChangeWrapper(index) {
    // Need to convert the index to the actual value in
    //   the dropdown
    onChange(list[index]);
    onBlur();
  }

  return (
    <Dropdown {...otherProps} list={list} onSelect={onChangeWrapper}></Dropdown>
  );
}
