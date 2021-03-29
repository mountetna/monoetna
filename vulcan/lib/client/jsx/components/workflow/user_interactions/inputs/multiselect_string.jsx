import React from 'react';

import ListInput from 'etna-js/components/inputs/list_input';
import DropdownInput from 'etna-js/components/inputs/dropdown_input';
import {nextUiStepsSelector} from '../../../../selectors/workflow';

export default function MultiselectStringInput({input, onChange}) {
  if (!input || !onChange) return null;

  var collator = new Intl.Collator(undefined, {
    numeric: true,
    sensitivity: 'base'
  });

  return (
    <ListInput
      placeholder='Select items from the list'
      className='link_text'
      values={
        input.default
          ? Array.isArray(input.default)
            ? input.default
            : [input.default]
          : []
      }
      itemInput={DropdownInput}
      list={(input.options && input.options.sort(collator.compare)) || []}
      onChange={(e) => {
        onChange(input.name, e);
      }}
    ></ListInput>
  );
}
