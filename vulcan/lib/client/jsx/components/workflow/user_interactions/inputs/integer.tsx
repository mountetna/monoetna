import React from 'react';
import {IntegerInput as EtnaIntegerInput} from 'etna-js/components/inputs/numeric_input';
import {InputBackendComponent} from './input_types';

const IntegerInput: InputBackendComponent = ({input, onChange}) => {
  if (!input || !onChange) return null;

  return (
    <EtnaIntegerInput
      defaultValue={input.value}
      onChange={(e: any) => {
        onChange(input.name, e);
      }}
    />
  );
};

export default IntegerInput;
