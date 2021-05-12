import React, {useState, useMemo, useEffect} from 'react';
import * as _ from 'lodash';

import {TYPE} from '../../../../api_types';
import {InputBackendComponent, InputSpecification} from './input_types';
import UserInput from './user_input';
import {inputValueNonEmpty} from '../../../../selectors/workflow_selectors';

const MultipleMultiselectStringAllInput: InputBackendComponent = ({
  input,
  onChange
}) => {
  interface Selection {
    cwlInputName: string;
    inputName: string;
    value: string[] | null;
  }

  const [selectedValues, setSelectedValues] = useState(
    {} as {[key: string]: {[key: string]: string[] | null}}
  );
  const [lastSelectedValues, setLastSelectedValues] = useState(
    {} as {[key: string]: {[key: string]: string[] | null}}
  );
  const [mockInputs, setMockInputs] = useState(
    {} as {[key: string]: InputSpecification[]}
  );
  const [userSelection, setUserSelection] = useState({} as Selection);

  const options = useMemo(() => {
    return input.data;
  }, [input.data]);

  function equalKeys(hash1: {}, hash2: {}): boolean {
    return _.isEqual(Object.keys(hash1), Object.keys(hash2));
  }

  function allSelectedValuesNonEmpty(hash: {}): boolean {
    return Object.values(hash).every((value) => inputValueNonEmpty(value));
  }

  const noEmptyStrings: boolean = useMemo(
    () =>
      Object.values(
        selectedValues
      ).every((nestedInput: {[key: string]: string[] | null}) =>
        Object.values(nestedInput).every((value) =>
          value?.every((val) => '' !== val)
        )
      ),
    [selectedValues]
  );

  // Helper methods to check if all inputs have values
  const allInputsPresent: boolean = useMemo(
    () => equalKeys(options || {}, selectedValues),
    [options, selectedValues]
  );

  const allNestedInputsPresent: boolean = useMemo(() => {
    if (!options || !allInputsPresent) return false;

    return Object.entries(options).every(([cwlInputName, nestedInput]) =>
      equalKeys(nestedInput, selectedValues[cwlInputName])
    );
  }, [options, selectedValues, allInputsPresent]);

  const allInputsPopulated: boolean = useMemo(() => {
    if (!options || !allNestedInputsPresent) return false;

    return Object.keys(options).every((cwlInputName) =>
      allSelectedValuesNonEmpty(selectedValues[cwlInputName])
    );
  }, [options, selectedValues, allNestedInputsPresent]);

  function selectedValuesChanged() {
    return !_.isEqual(selectedValues, lastSelectedValues);
  }

  useEffect(() => {
    let {cwlInputName, inputName, value} = userSelection;

    if (!cwlInputName || !inputName) return;

    setSelectedValues({
      ...selectedValues,
      [cwlInputName]: {
        ...selectedValues[cwlInputName],
        [inputName]: value
      }
    });
  }, [userSelection]);

  useEffect(() => {
    if (Object.keys(options || {}).length > 0 && selectedValuesChanged()) {
      if (allInputsPopulated) {
        if (noEmptyStrings) onChange(input.name, selectedValues);
      } else {
        onChange(input.name, null);
      }
    }

    setLastSelectedValues(selectedValues);
  }, [selectedValues]);

  useEffect(() => {
    if (options && input.default) {
      setSelectedValues(input.default);
    }
  }, [options]);

  useEffect(() => {
    if (options && selectedValues) {
      // We need to create "dummy" inputs with only a subset of the data,
      //      and the data key as the label.
      setMockInputs(
        Object.entries(options).reduce(
          (
            acc: {[key: string]: InputSpecification[]},
            [cwlInputName, nestedInputs]: [string, {[key: string]: string[]}]
          ): {[key: string]: InputSpecification[]} => {
            acc[cwlInputName] = Object.entries(nestedInputs).map(
              ([label, values]: [string, string[]]): InputSpecification => {
                return {
                  type: TYPE.MULTISELECT_STRING_ALL,
                  name: label,
                  label: label,
                  data: {[label]: values},
                  default:
                    selectedValues[cwlInputName] &&
                    selectedValues[cwlInputName][label]
                      ? selectedValues[cwlInputName][label]
                      : null
                };
              }
            );
            return acc;
          },
          {}
        )
      );
    }
  }, [options, selectedValues]);

  return (
    <div>
      {Object.keys(mockInputs).map((cwlInputName: string) => {
        let groupedInputs: InputSpecification[] = mockInputs[cwlInputName];
        return (
          <div key={cwlInputName}>
            {groupedInputs.map((mockInput) => {
              return (
                <UserInput
                  key={`${cwlInputName}-${mockInput.name}`}
                  onChange={(inputName, value) => {
                    setUserSelection({
                      cwlInputName,
                      inputName,
                      value
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
