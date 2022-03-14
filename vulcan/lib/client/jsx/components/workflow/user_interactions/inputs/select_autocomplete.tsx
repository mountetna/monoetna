import React, { useCallback, useEffect, useMemo, useState } from 'react';
import {DataEnvelope, WithInputParams} from './input_types';
import {maybeOfNullable, some, withDefault} from "../../../../selectors/maybe";
import {flattenStringOptions, StringOptions} from "./monoids";
import { useMemoized } from "../../../../selectors/workflow_selectors";
import { useSetsDefault } from './useSetsDefault';
import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';
import { debounce } from 'lodash';


export default function SelectAutocompleteInput(
  {data, label, minWidth, maxOptions=100, clearable=true, onChange, ...props}: WithInputParams<
  {label?: string, minWidth?: number, clearable?: boolean, maxOptions?:number}, string | null, StringOptions>) {
  /*
  Creates a searchable dropdown selection input box from concatenated values of the 'data' hash.
  Special Case: If any data key is "recommendation", a line of text will display the values of this recommendation to the user.
  */
  const [data_use, suggestion] = useMemoized(pullRecommendation, data)
  const options_in = useMemoized(flattenStringOptions, data_use);
  const value = useSetsDefault(null, props.value, onChange);
  const disp_label = useMemo( () => {return suggestion ? suggestion : label}, [suggestion, label]);
  const [options, setOptions] = useState(options_in as typeof options_in);
  const [loadingOptions, setLoadingOptions] = useState(false);
  const [helperText, setHelperText] = useState(undefined as string | undefined);
  const [inputState, setInputState] = React.useState(dispValue(value));
  const getOptionsDelayed = useCallback(
    debounce((text, options_in, callback) => {
      // console.log('inputState',inputState)
      // console.log('value',value)
      // console.log('options',options)
      setLoadingOptions(true); setOptions([]);
      getOptionsAsync(text, [...options_in]).then(callback);
    }, 200),
    [],
  );
  
  function determineHelperText(filteredOptions: typeof optionsLeaf) {
    if (filteredOptions.length>maxOptions) {
      setHelperText(filteredOptions.length + " options, only " + maxOptions + " shown")
    } else if (filteredOptions.length==0) {
      setHelperText("no matching options")
    } else {
      setHelperText(undefined)
    }
  }
  
  useEffect(() => {
    if (options_in.length>1000) {
      getOptionsDelayed(inputState, options_in, (filteredOptions: typeof options) => {
        console.log('calculating options - slow')
        setLoadingOptions(false)
        determineHelperText(filteredOptions)
        setOptions(filteredOptions.splice(0,maxOptions-1))
      });
    } else {
      console.log('calculating options - fast')
      const filteredOptions = filterOptions(inputState, options_in)
      determineHelperText(filteredOptions)
      setOptions(filteredOptions)
    }
  }, [inputState, getOptionsDelayed, options_in]);
  // console.log('inputState',inputState)
  // console.log('value',value)
  // console.log('options',options)
  return (
    <Autocomplete
      disableClearable={clearable}
      clearOnBlur={true}
      options={options}
      // disable filtering on client
      filterOptions={(x: typeof options) => x}
      loading={loadingOptions}
      value={value}
      onChange={(event: any, newInputState: string | null) => {
        // console.log('inputState',inputState)
        // console.log('value',value)
        onChange(maybeOfNullable(newInputState))
      }}
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState)
      }}
      ListboxProps={{style:{paddingTop: helperText ? 0 : 20}}}
      style={{minWidth: minWidth, paddingTop: disp_label ? 8 : 0}}
      renderInput={(params:any) => (
        <TextField 
          {...params}
          helperText={helperText}
          error={inputState!=dispValue(value)}
          label={disp_label}
          size="small"
          InputLabelProps={{ shrink: true }}/>
        )}
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
        const rec = Array.isArray(data['recommendation']) ? data['recommendation'].join(', ') : data['recommendation'];
        suggestion = "Recommendation: " + rec;
        delete data_use['recommendation']
      }
    return [data_use, suggestion]
    
  }
  return [data, undefined]
}

const getOptionsAsync = (query: string, options: string[]) => {
  return new Promise((resolve) => {
    setTimeout(() => {
      resolve(
        filterOptions(query, options),
      );
    }, 3000);
  });
};

function filterOptions(query: string, options: string[]) {
  return options.filter(
    (o) => {
      return (query==null) ? true : o.indexOf(query) > -1
    },
  )
}

function dispValue(value: string|null) {
  return (value==null) ? '' : value
}