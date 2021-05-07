import React, {useCallback, useState, useEffect} from 'react';

import MultiselectStringAllInput from './multiselect_string_all';
import {InputBackendComponent, InputSpecification} from './input_types';

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
    {} as {[key: string]: string | null}
  );
  const [mockInputs, setMockInputs] = useState([] as InputSpecification[]);

  const onChangeSingleInput = (inputName: string, value: string) => {
    setSelectedValues({...selectedValues, [inputName]: value});
  };

  const generateMockInputs = useCallback((): InputSpecification[] => {
    return Object.keys(input.data || {}).map(
      (key: string): InputSpecification => {
        return {
          type: 'multiselect_string_all',
          name: key,
          label: key,
          data: {[key]: input.data ? input.data[key] : []},
          default: null
        };
      }
    );
  }, [input]);

  useEffect(() => {
    setSelectedValues(
      Object.keys(input.data || {}).reduce(
        (acc: {[key: string]: null}, key: string) => {
          acc[key] = null;
          return acc;
        },
        {}
      )
    );
    setMockInputs(generateMockInputs());
  }, [input.data]);

  useEffect(() => {
    if (
      Object.values(selectedValues).filter((value) => null == value).length ===
      0
    ) {
      onChange(input.name, selectedValues);
    }
  }, [selectedValues]);

  return (
    <div>
      {mockInputs.map((mockInput) => {
        return (
          <MultiselectStringAllInput
            onChange={onChangeSingleInput}
            input={mockInput}
          />
        );
      })}
    </div>
  );
};

export default MultipleMultiselectStringAllInput;
