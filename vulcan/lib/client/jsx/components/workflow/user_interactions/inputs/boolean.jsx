import React from 'react';

export default function BooleanInput({input, onChange}) {
  if (!input || !onChange) return null;

  return (
    <input
      type='checkbox'
      className='text_box'
      onChange={(e) => {
        onChange(input.name, e);
      }}
      defaultChecked={input.default}
    />
  );
}
