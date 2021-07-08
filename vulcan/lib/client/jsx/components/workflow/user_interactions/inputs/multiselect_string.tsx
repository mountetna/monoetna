import React, {Dispatch, useCallback, useMemo} from 'react';

import ListInput from 'etna-js/components/inputs/list_input';
import DropdownAutocompleteInput from 'etna-js/components/inputs/dropdown_autocomplete_wrapper';
import {InputBackendComponent, InputSpecification, WithInputParams} from './input_types';
import {mapSome, some, withDefault} from "../../../../selectors/maybe";

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

export default function MultiselectStringInput({onChange, data, ...props}: WithInputParams<{onClear?: Dispatch<void>, onAll?: Dispatch<void>}>) {
  const defaultOnClear = useCallback(() => onChange(null), [onChange]);
  const options = useMemo(() => getAllOptions(data).sort(collator.compare), [data]);
  const defaultOnAll = useCallback(() => onChange(some(options)), [onChange, options]);
  const {onClear = defaultOnClear, onAll = defaultOnAll} = props;
  const value = withDefault(mapSome(props.value, inner => Array.isArray(inner) ? inner : [inner]), []);

  return (
    <ListInput
      placeholder='Select items from the list'
      className='link_text'
      values={value}
      itemInput={DropdownAutocompleteInput}
      list={options}
      onChange={(e: string[]) => onChange(some(e))}
      onAll={onAll}
      onClear={onClear}
      maxItems={25}
    />
  );
}