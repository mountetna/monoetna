import React, {useMemo} from 'react';

import ListInput from 'etna-js/components/inputs/list_input';
import DropdownAutocompleteInput from 'etna-js/components/inputs/dropdown_autocomplete_wrapper';
import {InputBackendComponent, InputSpecification} from './input_types';

export function getAllOptions(data: InputSpecification['data']): string[] {
  return Object.values(data || {}).reduce((acc, n) => {
    if (Array.isArray(n)) acc.push(...n);
    else if (typeof n === "string") acc.push(n);
    else acc.push(n + "");

    return acc;
  }, []);
}

const collator = new Intl.Collator(undefined, {
  numeric: true,
  sensitivity: 'base'
});

const MultiselectStringInput: InputBackendComponent = ({input, onChange, onClear, onAll}) => {
  const options = useMemo(() => getAllOptions(input.data).sort(collator.compare), [input.data]);

  return (
    <ListInput
      placeholder='Select items from the list'
      className='link_text'
      values={
        input.value
          ? Array.isArray(input.value)
            ? input.value
            : [input.value]
          : []
      }
      itemInput={DropdownAutocompleteInput}
      list={options}
      onChange={(e: any) => {
        onChange(input.name, e);
      }}
      onAll={onAll}
      onClear={onClear}
      maxItems={25}
    />
  );
}

export default MultiselectStringInput;