import React from 'react';
import SlowTextInput from 'etna-js/components/inputs/slow_text_input';
import {WithInputParams} from './input_types';
import {some} from "../../../../selectors/maybe";
import {useSetsDefault} from "./useSetsDefault";
import {selectDefaultString} from "./monoids";

export default function StringInput({onChange, label, data, ...props}: WithInputParams<{label?: string}, string, string>) {
  const value = useSetsDefault(selectDefaultString(data), props.value, onChange);

  const inner = <SlowTextInput
      followDefault
      defaultValue={value}
      onChange={(e: string) => {
        onChange(some(e));
      }}
    />

  if (label) {
  return <div>
    {label}
    {inner}
  </div>;
  }

  return inner;
}