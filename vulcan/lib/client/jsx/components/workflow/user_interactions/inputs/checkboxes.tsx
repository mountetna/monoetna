import React, {useCallback, useEffect, useMemo, useState} from 'react';
import {InputBackendComponent} from './input_types';
import {getAllOptions} from './multiselect_string';

function CheckboxInput({
  onChange,
  option,
  checked
}: {
  onChange: (v: any) => void;
  option: any;
  checked: boolean;
}) {
  return (
    <label className='checkbox-input-option'>
      <input
        type='checkbox'
        className='text_box'
        onChange={() => {
          onChange(option);
        }}
        checked={checked}
      />
      {option}
    </label>
  );
}

const CheckboxesInput: InputBackendComponent = ({input, onChange}) => {
  const options = useMemo(() => getAllOptions(input.data).sort(), [input.data]);

  useEffect(() => {
    // Set all options as checked, only if no input.value provided.
    // Otherwise set checked options to mirror input.value
    if (null == input.value && options.length > 0) {
      onChange(input.name, [...options]);
    }
  }, [options, input.value, input.name, onChange]);

  const handleClickOption = useCallback(
    (option: any) => {
      let copy = [...(input.value || [])];
      if (!input.value.includes(option)) {
        copy.push(option);
      } else {
        copy = input.value.filter((opt: string) => option !== opt);
      }
      onChange(input.name, copy);
    },
    [input.value, onChange, input.name]
  );

  if (!input || !onChange) return null;

  return (
    <div className='checkbox-input-wrapper'>
      {options.map((option, index) => {
        return (
          <CheckboxInput
            option={option}
            key={index}
            checked={input.value ? input.value.includes(option) : false}
            onChange={handleClickOption}
          />
        );
      })}
    </div>
  );
};

export default CheckboxesInput;
