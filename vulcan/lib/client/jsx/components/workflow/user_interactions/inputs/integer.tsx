import React, {useCallback} from 'react';
import {IntegerInput as EtnaIntegerInput} from 'etna-js/components/inputs/numeric_input';
import {WithInputParams} from './input_types';
import {some} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import {selectDefaultNumber} from "./monoids";

export default function IntegerInput({onChange, data, ...props}: WithInputParams<{}, number, number>) {
  const onNewInt = useCallback((f: number) => onChange(some(f)), [onChange])
  const value = useSetsDefault(selectDefaultNumber(data), props.value, onChange);

  return (
    <EtnaIntegerInput
      followDefault
      defaultValue={value}
      onChange={onNewInt}
    />
  );
};