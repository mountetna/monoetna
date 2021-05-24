// Input component that takes a simple object and
// shows values based on a selected key

import React, {useState, useEffect, useMemo, useContext, useCallback} from 'react';
import * as _ from 'lodash';

import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import CheckboxesInput from './checkboxes';
import MultiselectStringInput, {getAllOptions} from './multiselect_string';
import {InputBackendComponent} from "./input_types";

const handleClickOption = useCallback((option: any) => {
  let copy = [...selectedOptions];
  if (!selectedOptions.includes(option)) {
    copy.push(option);
  } else {
    copy = selectedOptions.filter((opt) => option !== opt);
  }
  setSelectedOptions(copy);
}, [selectedOptions]);


const SingleDropdownMulticheckbox: InputBackendComponent = ({input, onChange}) => {
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
        <MultiselectStringInput
          onAll={handleAllInputs}
          onClear={handleClearInputs}
          onChange={onChange}
          input={input} />
      </div>
      <div className='checkbox-input-wrapper'>
        {options.map((option, index) => {
          return (
            <CheckboxInput
              option={option}
              key={index}
              checked={selectedOptions.includes(option)}
              onChange={handleClickOption}
            />
          );
        })}
      </div>
    </div>
  );
}

export default SingleDropdownMulticheckbox;
