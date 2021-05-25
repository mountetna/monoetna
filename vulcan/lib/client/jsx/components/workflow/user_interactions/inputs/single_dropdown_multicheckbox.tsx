// Input component that takes a simple object and
// shows values based on a selected key

import React, {useState, useEffect, useMemo, useCallback} from 'react';

import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import CheckboxesInput from './checkboxes';
import {InputBackendComponent, InputSpecification} from './input_types';

const SingleDropdownMulticheckbox: InputBackendComponent = ({
  input,
  onChange
}) => {
  const [checkboxInputs, setCheckboxInputs] = useState(
    {} as {[key: string]: InputSpecification | null}
  );
  const [selected, setSelected] = useState(
    {} as {[key: string]: {[key: string]: string[]}}
  );
  const [dropdownValues, setDropdownValues] = useState(
    {} as {[key: string]: string}
  );
  const [initialized, setInitialized] = useState(false);
  const {data} = input;

  const dropdownOptions: {[key: string]: string[]} = useMemo(() => {
    console.log('here', data);
    return Object.keys(data || {}).reduce(
      (
        acc: {[key: string]: string[]},
        nextCwlInput: string
      ): {[key: string]: string[]} => {
        acc[nextCwlInput] = data ? Object.keys(data[nextCwlInput]) : [];
        return acc;
      },
      {}
    );
  }, [data]);

  useEffect(() => {
    console.log('dropdownOptions', dropdownOptions);
    // By default set the first option as selected for each CWL input.
    setDropdownValues(
      Object.entries(dropdownOptions).reduce(
        (
          acc: {[key: string]: string},
          [nextCwlInput, options]: [string, string[]]
        ) => {
          acc[nextCwlInput] = options[0];
          return acc;
        },
        {}
      )
    );

    if (input.default) {
      // Set the boxes to any previous selection
      setSelected(input.default);
    } else {
      // by default all boxes are selected
      setSelected(
        Object.entries(dropdownOptions).reduce(
          (
            acc: {[key: string]: {[key: string]: string[]}},
            [nextCwlInput, options]: [string, string[]]
          ): {[key: string]: {[key: string]: string[]}} => {
            acc[nextCwlInput] = {};
            options.forEach((option) => {
              acc[nextCwlInput][option] = data
                ? data[nextCwlInput][option]
                : [];
            });
            return acc;
          },
          {}
        )
      );
    }
  }, [data]);

  useEffect(() => {
    // when options change for the dropdown menu,
    //   set the checkboxes to the options
    //   for the selected dropdownOption.
    setCheckboxInputs(
      Object.entries(dropdownValues).reduce(
        (
          acc: {[key: string]: InputSpecification},
          [cwlInputName, dropdownValue]: [string, string]
        ): {[key: string]: InputSpecification} => {
          acc[cwlInputName] = {
            ...input,
            name: `${cwlInputName}-${dropdownValue}`,
            data: data ? data[cwlInputName][dropdownValue] : [],
            default:
              selected && selected[cwlInputName]
                ? selected[cwlInputName][dropdownValue]
                : []
          };
          return acc;
        },
        {}
      )
    );
  }, [dropdownValues, selected]);

  useEffect(() => {
    // Update state whenever a new selection comes in.
    onChange(input.name, selected);
  }, [selected]);

  const handleSelect = useCallback(
    (cwlInputName: string, value: string) => {
      // Set the checkbox input to `null` first, to clear
      //   out the component.
      setCheckboxInputs({
        ...checkboxInputs,
        [cwlInputName]: null
      });
      setDropdownValues({...dropdownValues, [cwlInputName]: value});
    },
    [dropdownOptions]
  );

  const handleClickOption = useCallback(
    (cwlInputName: string, options: string[]) => {
      setSelected({
        ...selected,
        [cwlInputName]: {
          ...selected[cwlInputName],
          [dropdownValues[cwlInputName]]: [...options]
        }
      });
    },
    [dropdownValues]
  );

  if (!input || !onChange) return null;
  console.log('dropdownOptions2', dropdownOptions);
  return (
    <div>
      {Object.keys(dropdownValues).map((cwlInputName: string) => {
        let checkboxInput = checkboxInputs[cwlInputName];
        if (!checkboxInput) return null;

        return (
          <React.Fragment>
            <div>
              <DropdownAutocomplete
                defaultValue={dropdownValues[cwlInputName]}
                onSelect={(value: string) => handleSelect(cwlInputName, value)}
                list={dropdownOptions[cwlInputName]}
              />
            </div>
            <div className='checkbox-input-wrapper'>
              <CheckboxesInput
                key={`${cwlInputName}-${dropdownValues[cwlInputName]}`}
                input={checkboxInput}
                onChange={(inputName: string, checked: string[]) =>
                  handleClickOption(cwlInputName, checked)
                }
              />
            </div>
          </React.Fragment>
        );
      })}
    </div>
  );
};

export default SingleDropdownMulticheckbox;
