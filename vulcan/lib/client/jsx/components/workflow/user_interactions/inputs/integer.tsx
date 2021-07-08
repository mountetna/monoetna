import React, {useCallback} from 'react';
import {IntegerInput as EtnaIntegerInput} from 'etna-js/components/inputs/numeric_input';
import {InputBackendComponent} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";

const IntegerInput: InputBackendComponent = ({value, onChange}) => {
  const onNewInt = useCallback((f: number) => onChange(some(f)), [onChange])

  return (
    <EtnaIntegerInput
      defaultValue={withDefault(value, 0)}
      onChange={onNewInt}
    />
  );
};

export default IntegerInput;
