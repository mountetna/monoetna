import React, {useMemo, useCallback, useState, useEffect} from 'react';

import {TYPE} from '../../../../api_types';
import {InputBackendComponent, InputSpecification} from './input_types';
import UserInput from './user_input';

export const flattenOptions = (
  data: {[key: string]: {[key: string]: string[]}} | null | undefined
) => {
  if (data) {
    return Object.values(data).reduce((a, b) => ({...a, ...b}), {});
  }

  return {};
};

const MultipleMultiselectStringAllInput: InputBackendComponent = ({
  input,
  onChange
}) => {
  const {data, name, value: inputValue} = input;

  const options: {[k: string]: string[]} = useMemo(() => flattenOptions(data), [
    data
  ]);

  const innerInputs: InputSpecification[] = useMemo(() => {
    const result: InputSpecification[] = [];

    Object.keys(options).forEach((label) => {
      const innerValues = options[label];
      result.push({
        type: TYPE.MULTISELECT_STRING_ALL,
        label,
        name: label,
        data: {[label]: innerValues},
        value: inputValue ? inputValue[label] : null
      });
    });

    return result;
  }, [options, inputValue]);

  const onSelectInnerInput = useCallback(
    ({label, value}: {label: string; value: string[]}) => {
      const update: {[k: string]: string[] | null} = {
        ...inputValue,
        [label]: value
      };

      onChange(name, update);
    },
    [inputValue, name, onChange]
  );

  return (
    <div>
      {innerInputs.map((innerInput: InputSpecification) => {
        return (
          <UserInput
            key={innerInput.name}
            onChange={(label, value) => {
              onSelectInnerInput({
                label,
                value
              });
            }}
            input={innerInput}
          />
        );
      })}
    </div>
  );
};

export default MultipleMultiselectStringAllInput;
