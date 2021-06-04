import React, {useState, useMemo, useEffect, useCallback, useRef} from 'react';
import * as _ from 'lodash';

import {InputBackendComponent} from './input_types';
import TextInput from 'etna-js/components/inputs/text_input.jsx';

function equalKeys(hash1: {}, hash2: {}): boolean {
  return _.isEqual(Object.keys(hash1), Object.keys(hash2));
}

const MultipleStringInput: InputBackendComponent = ({
  input,
  onChange
}) => {
  const {data, name} = input;

  const options: {[k: string]: string} = useMemo(() => {
    if (data) {
      return Object.values(data).reduce((a, b) => ({...a, ...b}), {});
    };

    return {};
  }, [data]);

  useEffect(() => {
    // Set all key's values as the initially given values, only if no input.value provided. (Initialization)
    // Otherwise set to mirror input.value
    if (null == input.value) {
      onChange(name, options);
    }
  }, [options, input.value, name, onChange]);

  const updateLabel = (newValue: string, key: string, prevLabels = input.value) => {
    prevLabels[key] = newValue;
    onChange(name, prevLabels);
    console.log(input.value);
  };

  if (null == input.value) return null;

  return (
    <div>
      {Object.entries(input.value).map(([k, val]) => {
        return (
          <TextInput
            header={k}
            value={val || ''}
            onChange={(newValue: string) => updateLabel(newValue, k)}
          />
        );
      })}
    </div>
  );

};

export default MultipleStringInput;
