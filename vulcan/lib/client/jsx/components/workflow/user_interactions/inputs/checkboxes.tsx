import React, {useCallback, useEffect, useMemo, useState} from 'react';
import {InputBackendComponent, WithInputParams} from './input_types';
import {getAllOptions} from './multiselect_string';
import {some, withDefault} from "../../../../selectors/maybe";
import BooleanInput from "./boolean";

export default function CheckboxesInput({data, onChange, ...props}: WithInputParams<{}>) {
  const options = useMemo(() => getAllOptions(data).sort(), [data]);
  const value = useMemo(() => withDefault(props.value, []), [props.value]);

  useEffect(() => {
    // Set all options as checked, only if no input.value provided.
    // Otherwise set checked options to mirror input.value
    if (!value && options.length > 0) {
      onChange(some([...options]));
    }
  }, [options, onChange, value]);

  const handleClickOption = useCallback(
    (option: any) => {
      const newValue = [...value];
      const existingIndex = newValue.indexOf(option);
      if (existingIndex >= 0) newValue.splice(existingIndex, 1);
      else newValue.push(option);

      onChange(some(newValue));
    },
    [value, onChange]
  );

  return (
    <div className='checkbox-input-wrapper'>
      {options.map((option, index) => {
        return (
          <BooleanInput
            label={option}
            data={null}
            key={index}
            value={some(value.includes(option))}
            onChange={handleClickOption}
          />
        );
      })}
    </div>
  );
};
