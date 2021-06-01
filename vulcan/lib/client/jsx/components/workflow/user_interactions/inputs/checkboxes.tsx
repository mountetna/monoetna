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
  // because we want to load the widget as all inputs are
  //   checked, we'll still use local state here.
  const [selectedOptions, setSelectedOptions] = useState(input.value || []);

  const options = useMemo(() => getAllOptions(input.data).sort(), [input.data]);

  useEffect(() => {
    // on mount, set all options as checked
    setSelectedOptions([...options]);
  }, []);

  useEffect(() => {
    if (null != input.value) setSelectedOptions([...input.value]);
  }, [input.value]);

  const handleClickOption = useCallback(
    (option: any) => {
      let copy = [...selectedOptions];
      if (!selectedOptions.includes(option)) {
        copy.push(option);
      } else {
        copy = selectedOptions.filter((opt: string) => option !== opt);
      }
      onChange(input.name, copy);
    },
    [selectedOptions]
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
