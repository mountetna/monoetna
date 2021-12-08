import React from 'react';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {WithInputParams} from './input_types';
import {maybeOfNullable, some, withDefault} from "../../../../selectors/maybe";
import {flattenStringOptions, StringOptions} from "./monoids";
import {useMemoized} from "../../../../selectors/workflow_selectors";
import { useSetsDefault } from './useSetsDefault';

export default function SelectAutocompleteInput({data, onChange, ...props}: WithInputParams<{}, string | null, StringOptions>) {
  /*
  Creates a searchable dropdown selection input box from concatenated values of the 'data' hash.
  Special Case: If any data key is "recommendation", a line of text will display the values of this recommendation to the user.
  */
  const [data_use, suggestion] = useMemoized(pullRecommendation, data)
  const options = useMemoized(flattenStringOptions, data_use);
  const value = useSetsDefault("", props.value, onChange);

  return (
    <div>
      <DropdownAutocomplete
      onSelect={(e: string | null) => {
        onChange(maybeOfNullable(e));
      }}
      list={options}
      value={value}
      maxItems={30}
    />
    {suggestion}
    </div>
  );
};

function pullRecommendation<T>(data: T): [T, string | null] {
  
  if (data != null) {
  
    let data_use = {...data};
    let suggestion = null;
      if (Object.keys(data).includes('recommendation')) {
        // Allow multi-recommendation?
        const rec = Array.isArray(data['recommendation']) ? data['recommendation'][0] : data['recommendation'];
        suggestion = "Recommendation: " + rec;
        delete data_use['recommendation']
      }
    return [data_use, suggestion]
    
  }
  return [data, null]
}
