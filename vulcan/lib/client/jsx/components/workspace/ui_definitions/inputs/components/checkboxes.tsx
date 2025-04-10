import React, {useCallback} from 'react';
import {WithInputParams} from '../../input_types';
import {some} from '../../../../../selectors/maybe';
import {flattenStringOptions, StringOptions} from '../../monoids';
import {useMemoized} from '../../../../../selectors/workflow_selectors';
import {useSetsDefault} from '../../useSetsDefault';
import { CheckboxPieceRct } from '../pieces/user_input_pieces';

export default function CheckboxesInput({data, defaultValue, onChange, ...props}: WithInputParams<{}, string[], StringOptions>) {
  const options: string[] = useMemoized(flattenStringOptions, data);
  const value: string[] = useSetsDefault(defaultValue || options, props.value, onChange, 'picked');

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
      {options.length>0 && options.map((option, index) => {
        return (
          <div key={index}>
            <CheckboxPieceRct
              name={`checkbox-${index}`}
              changeFxn={(v, k?) => handleClickOption(option)}
              value={value.includes(option)}
              label={option}
            />
          </div>
        );
      })}
    </div>
  );
};
