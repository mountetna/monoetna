import React, {useState, useMemo, useEffect, useCallback, useRef} from 'react';
import * as _ from 'lodash';

import {TYPE} from '../../../../api_types';
import {InputBackendComponent, InputSpecification} from './input_types';
import UserInput from './user_input';
import {inputValueNonEmpty} from '../../../../selectors/workflow_selectors';
import {useBufferedInputState} from "./buffered_input_state";

function equalKeys(hash1: {}, hash2: {}): boolean {
  return _.isEqual(Object.keys(hash1), Object.keys(hash2));
}

const MultipleMultiselectStringAllInput: InputBackendComponent = ({
  input,
  onChange
}) => {
  const {data, name} = input;

  const [selectedValues, setSelectedValues] = useBufferedInputState<{[k: string]: string[] | null}>(input, {});

  const options: {[k: string]: string[]} = useMemo(() => {
    if (data) {
      return Object.values(data).reduce((a, b) => ({...a, ...b}), {});
    }

    return {};
  }, [data]);

  const innerInputs: InputSpecification[] = useMemo(() => {
    const result: InputSpecification[] = [];

    Object.keys(options).forEach(label => {
      const innerValues = options[label];
      result.push({
        type: TYPE.MULTISELECT_STRING_ALL,
        label,
        name: label,
        data: {[label]: innerValues},
        default: selectedValues[label] || null,
      });
    });

    return result;
  }, [options, selectedValues]);

  const onSelectInnerInput = useCallback(({label, value}: {label: string, value: string[]}) => {
    const update: {[k: string]: string[] | null} = {...selectedValues, [label]: value};

    const allInputsPresent = !!options && equalKeys(options, selectedValues);
    const noEmptyStrings = Object.values(update).every(selected => selected?.every(val => '' !== val));
    const allNestedInputsPresent =  allInputsPresent && Object.entries(options).every(([label]) => label in update);
    const allInputsPopulated = allNestedInputsPresent && Object.values(update).every((nestedInputArray) => inputValueNonEmpty(nestedInputArray));

    setSelectedValues(update);

    console.log('update to', update);
    if (allInputsPopulated && noEmptyStrings) {
      console.log('onChange', name, update);
      onChange(name, update);
    } else if (!allInputsPopulated) {
      console.log('onChange', name, null);
      onChange(name, null);
    }
  }, [selectedValues, options, setSelectedValues, name, onChange]);

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
