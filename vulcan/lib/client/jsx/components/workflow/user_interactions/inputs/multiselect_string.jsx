import React from 'react';

import ListInput from 'etna-js/components/inputs/list_input';
import DropdownInput from 'etna-js/components/inputs/dropdown_input';

export default function MultiselectStringInput({input, onChange}) {
  if (!input || !onChange) return null;

  return (
    <ListInput
      placeholder='Select items from the list'
      className='link_text'
      values={input.default || []}
      itemInput={DropdownInput}
      list={input.options || []}
      onChange={(e) => {
        onChange(input.name, e);
      }}
    />
  );
}
