// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useState, useEffect, useCallback} from 'react';
import * as _ from 'lodash';

import {WithInputParams} from "./input_types";
import {some} from "../../../../selectors/maybe";
import {useMemoized} from "../../../../selectors/workflow_selectors";
import {joinNesting} from "./monoids";
import {useSetsDefault} from "./useSetsDefault";
import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';
import { debounce } from 'lodash';

function getPath(options: OptionSet, leaf: string): string[] {
  for (let [key, value] of Object.entries(options)) {
    if (value && typeof value === 'object') {
      // Look one step ahead for our leaf nodes
      if (Object.keys(value).includes(leaf) && null == value[leaf])
        return [key, leaf];
      const path = getPath(value, leaf)
      if (path.length > 0) return [key, ...getPath(value, leaf)];
    } else if (null == value && key == leaf) {
      return [key];
    }
  }

  return [];
}

function getOptions(
  desiredPath: string[],
  optionSet: OptionSet
): string[] | null {
  if (null == desiredPath || desiredPath.length === 0)
    return Object.keys(optionSet);
  // _.at always returns an array, so unwrap it.
  let results = _.at(optionSet, desiredPath.join('.'))[0];
  if (results) return Object.keys(results);
  return null;
}

function LeafOptions({
  options_in,
  value,
  depth,
  handleSelect,
  maxOptions=100,
}: {
  options_in: string[] | null;
  value: string | null;
  depth: number;
  handleSelect: (value: string | null, depth: number) => void;
  maxOptions?: number;
}) {
  if (!options_in) return null;
  console.log('options to '+depth+' level:', options_in.length)
  const [optionsLeaf, setOptionsLeaf] = useState(options_in as typeof options_in);
  const [loadingOptions, setLoadingOptions] = useState(false);
  const [helperText, setHelperText] = useState(undefined as string | undefined);
  const [inputState, setInputState] = React.useState(dispValue(value));
  const getOptionsDelayed = useCallback(
    debounce((text, options_in, callback) => {
      // console.log('inputState',inputState)
      // console.log('value',value)
      // console.log('options',options)
      setLoadingOptions(true); setOptionsLeaf([]);
      getOptionsAsync(text, [...options_in]).then(callback);
    }, 200),
    [],
  );
  
  useEffect(() => {
    if (options_in.length>1000) {
      getOptionsDelayed(inputState, options_in, (filteredOptions: typeof optionsLeaf) => {
        console.log('calculating options - slow')
        setLoadingOptions(false)
        setOptionsLeaf(filteredOptions.splice(0,maxOptions-1))
      });
    } else {
      console.log('calculating options - fast')
      setOptionsLeaf(filterOptions(inputState, options_in))
    }
    if (optionsLeaf.length>maxOptions) {
      setHelperText(optionsLeaf.length + " options, only " + maxOptions + " shown")
    } else if (optionsLeaf.length==0) {
      setHelperText("no matching options")
    } else {
      setHelperText(undefined)
    }
  }, [inputState, getOptionsDelayed, options_in]);

  return(
    <Autocomplete
      key={`${depth}-${options_in.slice(0, 5).join('-')}`}
      // disablePortal
      disableClearable
      clearOnBlur={false}
      value={value}
      onChange={(event:any, e: string | null) => 
        handleSelect(e, depth)}
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => 
        setInputState(newInputState)}
      options={optionsLeaf}
      // disable filtering on client
      filterOptions={(x: typeof optionsLeaf) => x}
      loading={loadingOptions}
      renderInput={(params:any) => (
        <TextField 
          {...params}
          helperText={helperText}
          error={inputState!=dispValue(value)}
          size="small"/>
        )}
    />
  )
}

type OptionSet = {[k: string]: null | OptionSet};

export default function NestedSelectAutocompleteInput({ label, data, onChange, ...props }: WithInputParams<{label?: string}, string|null, OptionSet>) {
  const value = useSetsDefault(null, props.value, onChange)
  const allOptions = useMemoized(joinNesting, data);
  const [path, setPath] = useState([] as string[]);
  
  useEffect(() => {
    if (value && value != null) {
      const updatedPath = getPath(allOptions, value);
      setPath(updatedPath);
    }
  }, [allOptions, value])

  const handleSelect = useCallback(
    (value: string | null, depth: number) => {
      // User has not selected something...perhaps
      //   still typing?
      if (null == value) return;

      // figure out where the next selection comes from, check if we're picking a leaf node.
      const updatedPath = path.slice(0, depth);
      updatedPath.push(value);
      setPath(updatedPath);

      if (getOptions(updatedPath, allOptions) == null) {
        // If we are updating a leaf
        onChange(some(value));
      } else {
        // Otherwise a leaf has not been selected.
        onChange(null);
      }
    },
    [allOptions, path, onChange]
  );

  let lab
  if (label) {
    lab = label;
  } else {
    lab = null;
  }
  console.log('value', value)
  return (
    <div>
      {lab}
      <div>
        {path.map((value, index) => {
              const options = getOptions(path.slice(0, index), allOptions);
              return (
                <Autocomplete
                  key={index}
                  // disablePortal
                  disableClearable
                  clearOnBlur={false}
                  onChange={(event:any, e: string | null) => 
                    handleSelect(e, index)}
                  options={options}
                  value={value}
                  // style={{minWidth: minWidth, paddingTop: disp_label ? 8 : 0}}
                  renderInput={(params:any) => (
                    <TextField 
                      {...params}
                      size="small"/>
                    )}
                />
              );
            })
        }
      </div>
      <div>
        <LeafOptions
          options_in={getOptions(path, allOptions)}
          value={value}
          depth={path.length}
          handleSelect={handleSelect}
        />
      </div>
    </div>
  );
};

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