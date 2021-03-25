// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useState, useEffect} from 'react';
import * as _ from 'lodash';

import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';

function LeafOptions({options, depth, handleSelect}) {
  if (!options) return null;
  return (
    <DropdownAutocomplete
      key={`${depth}-${options.slice(0, 5).join('-')}`}
      onSelect={(e) => {
        handleSelect(e, depth);
      }}
      list={options || []}
      defaultValue={null}
    ></DropdownAutocomplete>
  );
}

export default function NestedSelectAutocompleteInput({input, onChange}) {
  const [path, setPath] = useState([]);
  const [currentDepth, setCurrentDepth] = useState(-1);
  const [options, setOptions] = useState([]);
  const [initialized, setInitialized] = useState(false);
  const [originalOptions, setOriginalOptions] = useState(null);

  if (!input || !onChange) return null;

  function getPath(options, leaf) {
    for (let [key, value] of Object.entries(options)) {
      if (value && typeof value === 'object') {
        // Look one step ahead for our leaf nodes
        if (Object.keys(value).includes(leaf) && null == value[leaf])
          return [key, leaf];
        let path = getPath(value, leaf);
        if (path) return [key, ...path];
      } else if (null == value && key == leaf) {
        return [key];
      }
    }
  }

  useEffect(() => {
    if (!initialized && input.options.length > options.length) {
      // input.options is always an array, as part of
      //   the input component interface, so pull out the
      //   nested hash of options.
      let allOptions = input.options[0];
      setOriginalOptions({...allOptions});

      // If a default value is passed in, we need to
      //   find it and set it. Will have to assume
      //   that all leaves are unique, with only
      //   a singular path that leads to the leaf.
      if (input.default) {
        let path = getPath(allOptions, input.default);
        setPath(path);
        setOptions(getOptions(path, allOptions));
        setCurrentDepth(path.length - 1);
      } else {
        setOptions(Object.keys({...allOptions}));
      }
      setInitialized(true);
    }
  }, [input]);

  useEffect(() => {
    if (null == options) {
      // Return the leaf
      onChange(input.name, path[path.length - 1]);
    } else {
      // "unselect" this input
      onChange(input.name, null);
    }
  }, [path]);

  function getOptions(desiredPath, optionSet = originalOptions) {
    if (null == desiredPath) return Object.keys(optionSet);

    // _.at always returns an array, so unwrap it.
    let results = _.at(optionSet, desiredPath.join('.'))[0];

    if (results) return Object.keys(results);

    return null;
  }

  function handleSelect(value, depth) {
    // User has not selected something...perhaps
    //   still typing?
    if (null == value) return;

    let updatedPath;
    if (depth <= currentDepth) {
      updatedPath = path.slice(0, depth);
    } else {
      updatedPath = [...path];
    }

    updatedPath.push(value);
    setCurrentDepth(depth);
    setPath(updatedPath);
    setOptions(getOptions(updatedPath));
  }

  return (
    <div>
      <div>
        {path.length > 0
          ? path.map((value, index) => {
              return (
                <DropdownAutocomplete
                  key={index}
                  onSelect={(e) => {
                    handleSelect(e, index);
                  }}
                  list={getOptions(0 === index ? null : path.slice(0, index))}
                  defaultValue={value}
                ></DropdownAutocomplete>
              );
            })
          : null}
      </div>
      <div>
        <LeafOptions
          options={options}
          depth={currentDepth + 1}
          handleSelect={handleSelect}
        ></LeafOptions>
      </div>
    </div>
  );
}
