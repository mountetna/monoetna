import React, {useCallback, useState, useMemo, useEffect} from 'react';
import * as _ from 'lodash';

import {TYPE} from '../../../../api_types';
import {InputBackendComponent, InputSpecification} from './input_types';
import UserInput from './user_input';

const MultipleMultiselectStringAllInput: InputBackendComponent = ({
  input,
  onChange
}) => {
  // Here we do two things before generating the list of
  //   multiselects:
  // 1) Create "dummy" inputs with only a subset of the data,
  //      and the data key as the label.
  // 2) Create a "onChangeSingleInput" callback to only call
  //      onChange when all input keys have valid selections.
  const [selectedValues, setSelectedValues] = useState(
    {} as {[key: string]: {[key: string]: string[] | null}}
  );
  const [mockInputs, setMockInputs] = useState(
    {} as {[key: string]: InputSpecification[]}
  );
  const [lastInputValue, setLastInputValue] = useState(
    {} as {[key: string]: any} // really any should be string | string[] | null
  );

  const options = useMemo(() => {
    return input.data;
  }, [input.data]);

  // Helper methods to check if all inputs have values
  const allInputsPresent: boolean = useMemo(
    () => _.isEqual(Object.keys(options || {}), Object.keys(selectedValues)),
    [options, selectedValues]
  );

  const allNestedInputsPresent: boolean = useMemo(() => {
    if (!options || !allInputsPresent) return false;
    return Object.keys(options).every((cwlInputName) => {
      let nestedOption = options[cwlInputName];
      return _.isEqual(
        Object.keys(nestedOption),
        Object.keys(selectedValues[cwlInputName])
      );
    });
  }, [options, selectedValues, allInputsPresent]);

  const allInputsPopulated: boolean = useMemo(() => {
    if (!options || !allNestedInputsPresent) return false;
    return Object.keys(options).every((cwlInputName) => {
      return Object.values(selectedValues[cwlInputName]).every(
        (selectedValue) => null != selectedValue
      );
    });
  }, [options, selectedValues, allNestedInputsPresent]);

  useEffect(() => {
    let {cwlInputName, inputName, newValue} = lastInputValue;

    if (!cwlInputName || !inputName) return;

    setSelectedValues({
      ...selectedValues,
      [cwlInputName]: {
        ...selectedValues[cwlInputName],
        [inputName]: newValue
      }
    });
  }, [lastInputValue]);

  useEffect(() => {
    if (!options) return;

    setMockInputs(
      Object.keys(options).reduce(
        (
          acc: {[key: string]: InputSpecification[]},
          cwlInputName: string
        ): {[key: string]: InputSpecification[]} => {
          let subInputs: {[key: string]: string[]} = options[cwlInputName];

          acc[cwlInputName] = Object.keys(subInputs).map(
            (key: string): InputSpecification => {
              return {
                type: TYPE.MULTISELECT_STRING_ALL,
                name: key,
                label: key,
                data: {[key]: subInputs[key]},
                default:
                  selectedValues[cwlInputName] &&
                  selectedValues[cwlInputName][key]
                    ? selectedValues[cwlInputName][key]
                    : null
              };
            }
          );
          return acc;
        },
        {}
      )
    );
    
    if (allInputsPopulated) {
      onChange(input.name, selectedValues);
    } else {
      onChange(input.name, null);
    }
  }, [selectedValues]);

  useEffect(() => {
    if (options && input.default) {
      setSelectedValues(input.default);
    }
  }, [options]);

  return (
    <div>
      {Object.keys(mockInputs).map((cwlInputName: string) => {
        let groupedInputs: InputSpecification[] = mockInputs[cwlInputName];
        return (
          <div>
            {groupedInputs.map((mockInput) => {
              return (
                <UserInput
                  key={`${cwlInputName}-${mockInput.name}`}
                  onChange={(inputName, newValue) => {
                    setLastInputValue({
                      cwlInputName,
                      inputName,
                      newValue
                    });
                  }}
                  input={mockInput}
                />
              );
            })}
          </div>
        );
      })}
    </div>
  );
};

export default MultipleMultiselectStringAllInput;
