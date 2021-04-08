import React from 'react';
import {FloatInput as EtnaFloatInput} from 'etna-js/components/inputs/numeric_input';
import {InputBackendComponent} from "./input_types";

const FloatInput: InputBackendComponent = ({input, onChange}) => {
  if (!input || !onChange) return null;

  return (
    <EtnaFloatInput
      defaultValue={input.default}
      onChange={(e: any) => {
        onChange(input.name, e);
      }}
    />
  );
}

export default FloatInput;
