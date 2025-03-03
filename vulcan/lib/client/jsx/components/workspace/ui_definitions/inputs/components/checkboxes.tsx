import React, {useCallback} from 'react';
import {WithInputParams} from '../../input_types';
import {some} from '../../../../../selectors/maybe';
import BooleanInput from './boolean';
import {flattenStringOptions, StringOptions} from '../../monoids';
import {useMemoized} from '../../../../../selectors/workflow_selectors';
import {useSetsDefault} from '../../useSetsDefault';

export default function CheckboxesInput({data, defaultValue, onChange, ...props}: WithInputParams<{}, string[], StringOptions>) {
  const options = useMemoized(flattenStringOptions, data);
  const value = useSetsDefault(defaultValue || options, props.value, onChange, 'picked');

  const handleClickOption = useCallback(
    (option: string) => {
      const newValue = [...value];
      const existingIndex = newValue.indexOf(option);
      if (existingIndex >= 0) newValue.splice(existingIndex, 1);
      else newValue.push(option);
      onChange({picked: some(newValue)});
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
            labelPlacement='end'
            value={{value: some(value.includes(option))}}
            onChange={() => handleClickOption(option)}
          />
        );
      })}
    </div>
  );
};
