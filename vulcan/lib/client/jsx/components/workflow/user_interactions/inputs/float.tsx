import React, {useCallback} from 'react';
import {FloatInput as EtnaFloatInput} from 'etna-js/components/inputs/numeric_input';
import {WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";

export function FloatInput({value, onChange}: WithInputParams<{}, number>) {
  const onNewFloat = useCallback((f: number) => onChange(some(f)), [onChange])

  return (
    <EtnaFloatInput
      defaultValue={withDefault(value, 0)}
      onChange={onNewFloat}
    />
  );
}

export default FloatInput;
