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
  options,
  depth,
  handleSelect
}: {
  options: string[] | null;
  depth: number;
  handleSelect: (value: string | null, depth: number) => void;
}) {
  if (!options) return null;

  return(
    <Autocomplete
      key={`${depth}-${options.slice(0, 5).join('-')}`}
      disablePortal
      disableClearable
      onChange={(event:any, e: string | null) => 
        handleSelect(e, depth)}
      options={options}
      // style={{minWidth: minWidth, paddingTop: disp_label ? 8 : 0}}
      renderInput={(params:any) => (
        <TextField 
          {...params}
          variant="outlined"
          size="small"/>
        )}
    />
  )
}

type OptionSet = {[k: string]: null | OptionSet};

export default function NestedSelectAutocompleteInput({ label, data, onChange, ...props }: WithInputParams<{label?: string}, string, OptionSet>) {
  const value = useSetsDefault("", props.value, onChange)
  const allOptions = useMemoized(joinNesting, data);
  const [path, setPath] = useState([] as string[]);
  
  useEffect(() => {
    if (value && value != "") {
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

  return (
    <div>
      {lab}
      <div>
        {path.map((value, index) => {
              const options = getOptions(path.slice(0, index), allOptions);
              return (
                <Autocomplete
                  key={index}
                  disablePortal
                  disableClearable
                  onChange={(event:any, e: string | null) => 
                    handleSelect(e, index)}
                  options={options}
                  value={value}
                  // style={{minWidth: minWidth, paddingTop: disp_label ? 8 : 0}}
                  renderInput={(params:any) => (
                    <TextField 
                      {...params}
                      variant="outlined"
                      size="small"/>
                    )}
                />
              );
            })
        }
      </div>
      <div>
        <LeafOptions
          options={getOptions(path, allOptions)}
          depth={path.length}
          handleSelect={handleSelect}
        />
      </div>
    </div>
  );
};