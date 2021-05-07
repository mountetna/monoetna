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

  const options = useMemo(() => {
    return input.data;
  }, [input.data]);

  // Helper methods to check if all inputs have values
  const allInputsPresent: boolean = useMemo(
    () => _.isEqual(Object.keys(options || {}), Object.keys(selectedValues)),
    [options, selectedValues]
  );
  const allNestedInputsPresent: boolean = useMemo(() => {
    if (!options) return false;
    return (
      allInputsPresent &&
      Object.keys(options).every((cwlInputName) => {
        let nestedOption = options[cwlInputName];
        return _.isEqual(
          Object.keys(nestedOption),
          Object.keys(selectedValues[cwlInputName])
        );
      })
    );
  }, [options, selectedValues, allInputsPresent]);
  const allInputsPopulated: boolean = useMemo(() => {
    if (!options) return false;
    return (
      allNestedInputsPresent &&
      Object.keys(options).every((cwlInputName) => {
        return Object.values(selectedValues[cwlInputName]).every(
          (selectedValue) => null != selectedValue
        );
      })
    );
  }, [options, selectedValues, allNestedInputsPresent]);

  const onChangeSingleInput = useCallback(
    (cwlInputName: string, inputName: string, value: string[] | null) => {
      console.log('changing single input');
      console.log(
        'cwlInputName',
        cwlInputName,
        'inputName',
        inputName,
        'value',
        value
      );
      console.log('selectedValues', selectedValues);
      console.log({
        ...selectedValues,
        [cwlInputName]: {
          ...selectedValues[cwlInputName],
          [inputName]: value
        }
      });
      setSelectedValues({
        ...selectedValues,
        [cwlInputName]: {
          ...selectedValues[cwlInputName],
          [inputName]: value
        }
      });
    },
    [selectedValues]
  );

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
    console.log('allInputsPopulated', allInputsPopulated);
    if (allInputsPopulated) {
      onChange(input.name, selectedValues);
    } else {
      console.log('sending up the stack null');
      onChange(input.name, null);
    }
  }, [selectedValues]);

  useEffect(() => {
    console.log('options changed', selectedValues);
    if (options && input.default) {
      console.log('setting to default');
      setSelectedValues(input.default);
    }
  }, [options]);

  //   console.log('selectedValues', selectedValues);
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
                    onChangeSingleInput(cwlInputName, inputName, newValue);
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
