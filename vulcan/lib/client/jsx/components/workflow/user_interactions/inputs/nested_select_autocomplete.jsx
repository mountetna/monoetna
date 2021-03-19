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
      key={`${depth}-${options.join('-')}`}
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

  useEffect(() => {
    if (!initialized && input.options.length > options.length) {
      // input.options is always an array, as part of
      //   the input component interface, so pull out the
      //   nested hash of options.
      setOriginalOptions({...input.options[0]});
      setOptions(Object.keys({...input.options[0]}));
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
  }, [options]);

  function getOptions(desiredPath) {
    if (null == desiredPath) return Object.keys(originalOptions);

    // _.at always returns an array, so unwrap it.
    let results = _.at(originalOptions, desiredPath.join('.'))[0];

    if (results) return Object.keys(results);

    return null;
  }

  function handleSelect(value, depth) {
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
