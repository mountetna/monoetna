import React, {useState} from 'react';
import ListInput from 'etna-js/components/inputs/list_input';
import DropdownInput from 'etna-js/components/inputs/dropdown_input';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {
  IntegerInput,
  FloatInput
} from 'etna-js/components/inputs/numeric_input';
import SlowTextInput from 'etna-js/components/inputs/slow_text_input';

import {TYPE} from '../../../models/steps';
import InputHelp from './input_help';
import {wrapPaneItem} from '../../../utils/workflow';

export default function UserInput({input, onChange}) {
  const INPUTS = {
    default: (
      <SlowTextInput
        defaultValue={input.default}
        onChange={(e) => {
          onChange(input.name, e);
        }}
      ></SlowTextInput>
    ),
    [TYPE.INTEGER]: (
      <IntegerInput
        defaultValue={input.default}
        onChange={(e) => {
          onChange(input.name, e);
        }}
      ></IntegerInput>
    ),
    [TYPE.FLOAT]: (
      <FloatInput
        defaultValue={input.default}
        onChange={(e) => {
          onChange(input.name, e);
        }}
      ></FloatInput>
    ),
    [TYPE.BOOL]: (
      <input
        type='checkbox'
        className='text_box'
        onChange={(e) => {
          onChange(input.name, e);
        }}
        defaultChecked={input.default}
      />
    ),
    [TYPE.MULTISELECT_STRING]: (
      <ListInput
        placeholder='Select items from the list'
        className='link_text'
        values={input.default || []}
        itemInput={DropdownInput}
        list={input.options || []}
        onChange={(e) => {
          onChange(input.name, e);
        }}
      />
    ),
    [TYPE.SELECT_AUTOCOMPLETE]: (
      <DropdownAutocomplete
        onSelect={(e) => {
          onChange(input.name, e);
        }}
        list={input.options || []}
        defaultValue={input.default || null}
      ></DropdownAutocomplete>
    )
  };

  let inputType = input.type || 'default';

  let Component = <InputHelp input={input}>{INPUTS[inputType]}</InputHelp>;

  return wrapPaneItem(
    {
      name: input.name,
      value: Component
    },
    input.name
  );
}
