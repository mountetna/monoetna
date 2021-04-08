// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useState, useEffect, useMemo, useContext, useCallback} from 'react';
import * as _ from 'lodash';

import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {InputBackendComponent} from "./input_types";

function getPath(options: OptionSet, leaf: string): string[] {
  for (let [key, value] of Object.entries(options)) {
    if (value && typeof value === 'object') {
      // Look one step ahead for our leaf nodes
      if (Object.keys(value).includes(leaf) && null == value[leaf])
        return [key, leaf];
      let path = getPath(value, leaf);
      if (path.length > 0) return [key, ...path];
    } else if (null == value && key == leaf) {
      return [key];
    }
  }

  return [];
}

function getOptions(desiredPath: string[], optionSet: OptionSet): string[] | null {
  if (null == desiredPath || desiredPath.length === 0) return Object.keys(optionSet);
  // _.at always returns an array, so unwrap it.
  let results = _.at(optionSet, desiredPath.join('.'))[0];
  if (results) return Object.keys(results);
  return null;
}

function LeafOptions({options, depth, handleSelect}: {options: string[] | null, depth: number, handleSelect: (value: string | null, depth: number) => void}) {
  if (!options) return null;

  return (
    <DropdownAutocomplete
      key={`${depth}-${options.slice(0, 5).join('-')}`}
      onSelect={(e: any) => {
        handleSelect(e, depth);
      }}
      list={options}
      defaultValue={null}
    />
  );
}

type OptionSet = {[k: string]: null | OptionSet};

const NestedSelectAutocompleteInput: InputBackendComponent = ({input, onChange}) => {
  const [path, setPath] = useState([] as string[]);
  const {data} = input;

  if (!input || !onChange) return null;

  const allOptions: OptionSet = useMemo(() => {
    return Object.keys(data || {}).reduce((dataObj, k) => {
      if (!data) return dataObj;
      const next = data[k];
      return {
          ...dataObj,
          ...(typeof next === "object" && next != null && !Array.isArray(next) ? next : {})
      };
    }, {} as OptionSet);
  }, [data]);

  useEffect(() => {
    if (input.default) {
      setPath(getPath(allOptions, input.default));
    } else {
      setPath([])
    }
  }, [input, allOptions]);

  const handleSelect = useCallback((value: string | null, depth: number) => {
    // User has not selected something...perhaps
    //   still typing?
    if (null == value) return;

    const updatedPath = path.slice(0, depth);
    updatedPath.push(value);
    setPath(updatedPath);

    if (getOptions(updatedPath, allOptions) == null) {
      // If we are updating a leaf
      onChange(input.name, value);
    } else {
      // Otherwise a leaf has not been selected.
      onChange(input.name, null);
    }
  }, [input, allOptions, path, setPath, onChange]);

  return (
    <div>
      <div>
        {path.length > 0
          ? path.map((value, index) => {
            const options = getOptions(path.slice(0, index), allOptions);
              return (
                <DropdownAutocomplete
                  key={index}
                  onSelect={(e: string | null) => {
                    handleSelect(e, index);
                  }}
                  list={options}
                  defaultValue={value}
                />
              );
            })
          : null}
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
}

export default NestedSelectAutocompleteInput;
