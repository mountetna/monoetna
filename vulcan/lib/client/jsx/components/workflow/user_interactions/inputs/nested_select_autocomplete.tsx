// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useState, useEffect, useCallback} from 'react';
import * as _ from 'lodash';
import {WithInputParams} from './input_types';
import {some} from '../../../../selectors/maybe';
import {useMemoized} from '../../../../selectors/workflow_selectors';
import {joinNesting} from './monoids';
import {useSetsDefault} from './useSetsDefault';
import SelectAutocompleteInput from './select_autocomplete';

function getPath(options: OptionSet, leaf: string): string[] {
  for (let [key, value] of Object.entries(options)) {
    if (value && typeof value === 'object') {
      // Look one step ahead for our leaf nodes
      if (Object.keys(value).includes(leaf) && null == value[leaf])
        {return [key, leaf];}
      const path = getPath(value, leaf);
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
    {return Object.keys(optionSet);}
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
  label=undefined
}: {
  options_in: string[] | null;
  value: string | null;
  depth: number;
  handleSelect: (value: string | null, depth: number) => void;
  maxOptions?: number;
  label?: string | undefined
}) {
  if (!options_in) return null;
  if (value!=null && !options_in.includes(value)) {
    // non-null starting value, but path has not yet been properly determined?
    // console.log('skipping LeafOptions due to value mismatch')
    return null;
  }
  return (
    <SelectAutocompleteInput
      key={`${depth}-${options_in.slice(0, 5).join('-')}`}
      onChange={(v) => {}}
      onChangeOverride={
        (event:any, e: string|null) => handleSelect(e, depth)
      }
      minWidth={200}
      value={some(value)}
      data={{a: options_in}}
      maxOptions={maxOptions}
      label={label}
    />
  );
}

type OptionSet = {[k: string]: null | OptionSet};

export default function NestedSelectAutocompleteInput({ label, data, onChange, ...props }: WithInputParams<{label?: string}, string|null, OptionSet>) {
  const value = useSetsDefault(null, props.value, onChange);
  const allOptions = useMemoized(joinNesting, data);
  const [path, setPath] = useState([] as string[]);
  
  useEffect(() => {
    if (value!=null) {
      const updatedPath = getPath(allOptions, value);
      setPath(updatedPath);
    }
  }, [allOptions, value]);

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

  // console.log({value})
  // console.log({path})
  return (
    <div>
      {/* {lab} */}
      <div>
        {path.map((value, index) => {
              const options = getOptions(path.slice(0, index), allOptions);
              return (
                (options==null) ? null : // skip here if leaf
                <SelectAutocompleteInput
                  key={index}
                  onChange={(v) => {}}
                  onChangeOverride={
                    (event:any, e: string|null) => handleSelect(e, index)
                  }
                  minWidth={200}
                  value={some(value)}
                  data={{a: options}}
                  label={index==0 ? label : undefined}
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
          label={path.length == 0 ? label : undefined}
          handleSelect={handleSelect}
        />
      </div>
    </div>
  );
};
