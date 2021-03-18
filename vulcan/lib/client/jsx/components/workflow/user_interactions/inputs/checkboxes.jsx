import React, {useEffect, useState} from 'react';

function CheckboxInput({onChange, option}) {
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

export default function CheckboxesInput({input, onChange}) {
  const [selectedOptions, setSelectedOptions] = useState([]);
  const [initialized, setInitialized] = useState(false);

  if (!input || !onChange) return null;

  useEffect(() => {
    if (input.options.length > selectedOptions.length && !initialized) {
      setSelectedOptions([...input.options]);
      setInitialized(true);
    }
  }, [input]);

  useEffect(() => {
    onChange(input.name, selectedOptions);
  }, [selectedOptions]);

  function handleClickOption(option) {
    let copy = [...selectedOptions];
    if (!selectedOptions.includes(option)) {
      copy.push(option);
    } else {
      copy = selectedOptions.filter((opt) => option !== opt);
    }
    console.log(copy);
    setSelectedOptions(copy);
  }

  return (
    <div className='checkbox-input-wrapper'>
      {input.options.map((option, index) => {
        return (
          <CheckboxInput
            option={option}
            key={index}
            onChange={handleClickOption}
          ></CheckboxInput>
        );
      })}
    </div>
  );
}
