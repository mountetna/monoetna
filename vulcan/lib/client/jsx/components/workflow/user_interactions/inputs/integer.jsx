import React from 'react';
import {IntegerInput as EtnaIntegerInput} from 'etna-js/components/inputs/numeric_input';

export default function IntegerInput({input, onChange}) {
  if (!input || !onChange) return null;

  return (
    <EtnaIntegerInput
      defaultValue={input.default}
      onChange={(e) => {
        onChange(input.name, e);
      }}
    ></EtnaIntegerInput>
  );
}
