import React, {useCallback, useEffect, useMemo, useState} from 'react';
import {InputBackendComponent} from './input_types';
import {getAllOptions} from './multiselect_string';

import {useBufferedInputState} from './buffered_input_state';

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
  const [selectedOptions, setSelectedOptions] = useBufferedInputState<string[]>(
    input,
    []
  );

  const options = useMemo(() => getAllOptions(input.data).sort(), [input.data]);

  const handleClickOption = useCallback(
    (option: any) => {
      let copy = [...selectedOptions];
      if (!selectedOptions.includes(option)) {
        copy.push(option);
      } else {
        copy = selectedOptions.filter((opt) => option !== opt);
      }
      setSelectedOptions(copy);
      onChange(input.name, copy);
    },
    [selectedOptions, onChange]
  );

  if (!input || !onChange) return null;

  return (
    <div className='checkbox-input-wrapper'>
      {options.map((option, index) => {
        return (
          <CheckboxInput
            option={option}
            key={index}
            checked={selectedOptions.includes(option)}
            onChange={handleClickOption}
          />
        );
      })}
    </div>
  );
};

export default CheckboxesInput;
