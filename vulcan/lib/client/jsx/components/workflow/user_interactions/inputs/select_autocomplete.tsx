import React from 'react';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {DataEnvelope, WithInputParams} from './input_types';
import {maybeOfNullable, some, withDefault} from "../../../../selectors/maybe";
import {flattenStringOptions, StringOptions} from "./monoids";
import {useMemoized} from "../../../../selectors/workflow_selectors";
import { useSetsDefault } from './useSetsDefault';
import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';


export default function SelectAutocompleteInput({data, onChange, ...props}: WithInputParams<{}, string | null, StringOptions>) {
  /*
  Creates a searchable dropdown selection input box from concatenated values of the 'data' hash.
  Special Case: If any data key is "recommendation", a line of text will display the values of this recommendation to the user.
  */
  const [data_use, suggestion] = useMemoized(pullRecommendation, data)
  const options = useMemoized(flattenStringOptions, data_use);
  const value = useSetsDefault(null, props.value, onChange);

  return (
    <Autocomplete
      disablePortal
      value={value}
      onChange={(event:any, e: string | null) => 
        onChange(maybeOfNullable(e))}
      options={options}
      renderInput={(params:any) => <TextField {...params} label={suggestion} variant="outlined"/>}
    />
  );
};

function pullRecommendation<T extends DataEnvelope<any>>(data: T | null | undefined): [T | null | undefined, string | undefined] {
  // Make the recommendation part of the autocomplete into the 'label' of renderInput
  if (data != null) {
  
    let data_use = {...data};
    let suggestion = undefined;
      if (Object.keys(data).includes('recommendation')) {
        // Allow multi-recommendation?
        const rec = Array.isArray(data['recommendation']) ? data['recommendation'][0] : data['recommendation'];
        suggestion = "Recommendation: " + rec;
        delete data_use['recommendation']
      }
    return [data_use, suggestion]
    
  }
  return [data, undefined]
}
