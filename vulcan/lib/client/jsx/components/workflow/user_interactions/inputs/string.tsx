import React from 'react';
import SlowTextInput from 'etna-js/components/inputs/slow_text_input';
import {WithInputParams} from './input_types';
import {some, withDefault} from "../../../../selectors/maybe";

export default function StringInput({onChange, value}: WithInputParams<{}>) {
  return (
    <SlowTextInput
      defaultValue={withDefault(value, "")}
      onChange={(e: any) => {
        onChange(some(e));
      }}
    />
  );
}