import React, {useCallback} from 'react';
import {FloatInput as EtnaFloatInput} from 'etna-js/components/inputs/numeric_input';
import {WithInputParams} from './input_types';
import {some} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import {selectDefaultNumber} from "./monoids";

export function FloatInput({onChange, data, ...props}: WithInputParams<{}, number, number>) {
  const onNewFloat = useCallback((f: number) => onChange(some(f)), [onChange])
  const value = useSetsDefault(selectDefaultNumber(data), props.value, onChange);

  return (
    <EtnaFloatInput
      key={value}
      defaultValue={value}
      onChange={onNewFloat}
    />
  );
}

export default FloatInput;
