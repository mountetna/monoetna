import React from 'react';
import SlowTextInput from 'etna-js/components/inputs/slow_text_input';
import {InputBackendComponent} from "./input_types";

const StringInput: InputBackendComponent = ({input, onChange}) => {
  if (!input || !onChange) return null;

  return (
    <SlowTextInput
      defaultValue={input.default}
      onChange={(e: any) => {
        onChange(input.name, e);
      }}
    />
  );
}

export default StringInput;
