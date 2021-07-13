import React, {useCallback} from 'react';
import {IntegerInput as EtnaIntegerInput} from 'etna-js/components/inputs/numeric_input';
import {InputBackendComponent, WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";

export default function IntegerInput({value, onChange}: WithInputParams<{}, number>) {
  const onNewInt = useCallback((f: number) => onChange(some(f)), [onChange])

  return (
    <EtnaIntegerInput
      defaultValue={withDefault(value, 0)}
      onChange={onNewInt}
    />
  );
};