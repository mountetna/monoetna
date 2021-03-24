import React from 'react';
import SlowTextInput from 'etna-js/components/inputs/slow_text_input';

export default function StringInput({input, onChange}) {
  if (!input || !onChange) return null;

  return (
    <SlowTextInput
      defaultValue={input.default}
      onChange={(e) => {
        onChange(input.name, e);
      }}
    ></SlowTextInput>
  );
}
