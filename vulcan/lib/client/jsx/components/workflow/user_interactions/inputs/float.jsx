import React from 'react';
import {FloatInput as EtnaFloatInput} from 'etna-js/components/inputs/numeric_input';

export default function FloatInput({input, onChange}) {
  if (!input || !onChange) return null;

  return (
    <EtnaFloatInput
      defaultValue={input.default}
      onChange={(e) => {
        onChange(input.name, e);
      }}
    ></EtnaFloatInput>
  );
}
