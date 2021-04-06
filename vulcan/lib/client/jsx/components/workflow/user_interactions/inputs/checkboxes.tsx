import React, {useCallback, useEffect, useMemo, useState} from 'react';
import {InputBackendComponent} from "./input_types";
import {getAllOptions} from "./multiselect_string";

function CheckboxInput({onChange, option}: {onChange: (v: any) => void, option: any}) {
  return (
    <label className='checkbox-input-option'>
      <input
        type='checkbox'
        className='text_box'
        onChange={() => {
          onChange(option);
        }}
        defaultChecked={true}
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

    if (options.length > selectedOptions.length && !initialized) {
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
            onChange={handleClickOption}
          />
        );
      })}
    </div>
  );
}

export default CheckboxesInput;
