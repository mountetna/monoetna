import React, {useCallback, useEffect, useMemo, useState} from 'react';
import {InputBackendComponent} from "./input_types";
import {getAllOptions} from "./multiselect_string";

function CheckboxInput({onChange, option, checked}: {onChange: (v: any) => void, option: any, checked: boolean}) {
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
  const [selectedOptions, setSelectedOptions] = useState([] as any[]);
  const [initialized, setInitialized] = useState(false);
  if (!input || !onChange) return null;

  const options = useMemo(() => getAllOptions(input.data).sort(), [input.data]);

  useEffect(() => {
    // Setting any previously selected inputs (from storage or
    //   user interactions) takes precedence over setting
    //   all options as checked.
    if (input.default && input.default !== [] && !initialized) {
      setSelectedOptions([...input.default]);
      setInitialized(true);
    } else if (options.length > selectedOptions.length && !initialized) {
      setSelectedOptions([...options]);
      setInitialized(true);
    }
  }, [options]);

  useEffect(() => {
    onChange(input.name, selectedOptions);
  }, [selectedOptions]);

  const handleClickOption = useCallback((option: any) => {
    let copy = [...selectedOptions];
    if (!selectedOptions.includes(option)) {
      copy.push(option);
    } else {
      copy = selectedOptions.filter((opt) => option !== opt);
    }
    setSelectedOptions(copy);
  }, [selectedOptions]);

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
}

export default CheckboxesInput;
