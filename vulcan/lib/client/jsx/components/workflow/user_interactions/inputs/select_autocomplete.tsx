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
  {data, label, minWidth, maxOptions=100, disableClearable=true, onChangeOverride, onChange, ...props}: WithInputParams<
  {label?: string, minWidth?: number, disableClearable?: boolean, maxOptions?:number, onChangeOverride?: (event: any, e: string|null) => void}, string | null, StringOptions>) {
  /*
  Creates a searchable dropdown selection input box from concatenated values of the 'data' hash.
  Special Case: If any data key is "recommendation", a line of text will display the values of this recommendation to the user.
  */
  const [data_use, suggestion] = useMemoized(pullRecommendation, data)
  if (data_use == {}) return null
  const options_in = useMemoized(flattenStringOptions, data_use);
  const value = useSetsDefault(null, props.value, onChange);
  const disp_label = useMemo( () => {return suggestion ? suggestion : label}, [suggestion, label]);
  const [loadingOptions, setLoadingOptions] = useState(false);
  const [helperText, setHelperText] = useState(undefined as string | undefined);
  const [inputState, setInputState] = React.useState(dispValue(value));
  const [options, setOptionsFull] = useState({
    filtered: options_in,
    display: [...options_in].splice(0,maxOptions)
    });
  function setOptions(filteredOptions: string[]) {
    setOptionsFull({
      filtered: filteredOptions,
      display: [...filteredOptions].splice(0,maxOptions)})
  }
  const getOptionsDelayed = useCallback(
    debounce((text, options_in, callback) => {
      getOptionsAsync(text, [...options_in]).then(callback);
    }, 200),
    [],
  );
  
  useEffect(() => {
    // Shown from all options when user has made their selection (current text = current value)
    const query = (inputState!=value) ? inputState : ''
    
    if (options_in.length>1000) {
      setLoadingOptions(true);
      // console.log('calculating options - slow')
      getOptionsDelayed(query, options_in, (filteredOptions: string[]) => {
        // console.log('calculating options - slow')
        setLoadingOptions(false)
        setOptions(filteredOptions)
      });
    } else {
      // console.log('calculating options - fast')
      setOptions(filterOptions(query, options_in))
    }
  }, [inputState, getOptionsDelayed, options_in, value]);
  
  useEffect(() => {
    determineHelperText(setHelperText, options.filtered, options_in, loadingOptions, maxOptions, inputState, value)
  }, [options, options_in, inputState, value, maxOptions, loadingOptions, setHelperText])
  
  const onChangeAction = useCallback((event: any, e: string|null) => {
    (onChangeOverride) ? onChangeOverride(event, e) : onChange(maybeOfNullable(e))
  }, [onChangeOverride, onChange])
  
  return (
    <Autocomplete
      disableClearable={disableClearable}
      clearOnBlur={true}
      options={options_in}
      filterOptions={(x: string[]) => options.display}
      loading={loadingOptions}
      value={value}
      onChange={onChangeAction}
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState)
      }}
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

function dispValue(value: string|null) {
  return (value==null) ? '' : value
}

const getOptionsAsync = (query: string, opts: string[]) => {
  return new Promise((resolve) => {
    setTimeout(() => {
      resolve(
        filterOptions(query, opts),
      );
    }, 3000);
  });
};

function filterOptions(query: string, opts: string[]) {
  return opts.filter(
    (o) => {
      return (query==null) ? true : o.indexOf(query) > -1
    },
  )
}

function determineHelperText(
  setHelperText: React.Dispatch<React.SetStateAction<string | undefined>>,
  filteredOptions: string[],
  allOptions: string[],
  loadingOptions: boolean,
  maxOptions: number,
  typed: string,
  saved: string | null)
{
  
  if (loadingOptions) {
    return
  } else if (saved == null || typed!=saved) {
    if (filteredOptions.length>maxOptions) {
      setHelperText(allOptions.length + " options, only " + maxOptions + " shown")
    } else if (filteredOptions.length==0) {
      setHelperText("no matching options")
    } else {
      setHelperText(undefined)
    }
  } else {
    setHelperText(undefined)
  }
}