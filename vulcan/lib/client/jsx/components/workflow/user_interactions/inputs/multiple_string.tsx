import React, {useMemo, useCallback} from 'react';

import {WithInputParams} from './input_types';
import TextInput from 'etna-js/components/inputs/text_input';
import {withDefault} from "../../../../selectors/maybe";

export default function MultipleStringInput({ onChange, data, ...props }: WithInputParams<{}>) {
  const defaults: {[k: string]: string} = useMemo(() => {
    if (data) {
      return Object.values(data).reduce((a, b) => ({...a, ...b}), {});
    }

    return {};
  }, [data]);

  const value = withDefault(props.value, defaults);

  const updateLabel = useCallback((v: string, key: string) => {
    const newValue = {...value, [key]: v};
    onChange(newValue);
  }, [onChange, value])

  return (
    <div>
      {Object.entries(value).map(([k, val]) => {
        return (
          <TextInput
            key={k}
            header={k}
            value={val || ''}
            onChange={(newValue: string) => updateLabel(newValue, k)}
          />
        );
      })}
    </div>
  );

};