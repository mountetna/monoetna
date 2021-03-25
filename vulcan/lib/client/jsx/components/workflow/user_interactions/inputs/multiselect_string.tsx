import React from 'react';

import ListInput from 'etna-js/components/inputs/list_input';
import DropdownInput from 'etna-js/components/inputs/dropdown_input';
import {InputBackendComponent} from "./input_types";

const MultiselectStringInput: InputBackendComponent = ({input, onChange}) => {
  if (!input || !onChange) return null;

  const options: any[] = Object.values(input.data || {}).reduce((acc, n) => {
     if (Array.isArray(n)) return acc.concat(n);
     return acc.concat([n]);
  }, [input.data]);

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
      list={options}
      onChange={(e: any) => {
        onChange(input.name, e);
      }}
    />
  );
}

export default MultiselectStringInput;