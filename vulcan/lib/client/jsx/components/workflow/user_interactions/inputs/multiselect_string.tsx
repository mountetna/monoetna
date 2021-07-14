import React, {Dispatch, useCallback, useMemo} from 'react';

import ListInput from 'etna-js/components/inputs/list_input';
import DropdownAutocompleteInput from 'etna-js/components/inputs/dropdown_autocomplete_wrapper';
import {DataEnvelope, InputBackendComponent, BoundInputSpecification, WithInputParams} from './input_types';
import {mapSome, some, withDefault} from "../../../../selectors/maybe";
import {flattenStringOptions, StringOptions} from "./monoids";
import {useMemoized} from "../../../../selectors/workflow_selectors";

export default function MultiselectStringInput({onChange, data, ...props}: WithInputParams<
  {onClear?: Dispatch<void>, onAll?: Dispatch<void>}, string[], StringOptions>) {
  const defaultOnClear = useCallback(() => onChange(null), [onChange]);
  const options = useMemoized(flattenStringOptions, data);
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