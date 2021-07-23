import React from 'react';
import SlowTextInput from 'etna-js/components/inputs/slow_text_input';
import {WithInputParams} from './input_types';
import {some} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import {selectDefaultString} from "./monoids";

export default function StringInput({onChange, data, ...props}: WithInputParams<{}, string, string>) {
  const value = useSetsDefault(selectDefaultString(data), props.value, onChange);

  return (
    <SlowTextInput
      key={value}
      defaultValue={value}
      onChange={(e: string) => {
        onChange(some(e));
      }}
    />
  );
}