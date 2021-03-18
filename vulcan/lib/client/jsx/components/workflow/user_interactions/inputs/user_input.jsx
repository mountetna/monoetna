import React from 'react';

import {TYPE} from '../../../../models/steps';
import InputHelp from './input_help';
import BooleanInput from './boolean';
import FloatInput from './float';
import IntegerInput from './integer';
import MultiselectStringInput from './multiselect_string';
import SelectAutocompleteInput from './select_autocomplete';
import StringInput from './string';

export default function UserInput({input, onChange}) {
  const INPUTS = {
    default: StringInput,
    [TYPE.INTEGER]: IntegerInput,
    [TYPE.FLOAT]: FloatInput,
    [TYPE.BOOL]: BooleanInput,
    [TYPE.MULTISELECT_STRING]: MultiselectStringInput,
    [TYPE.SELECT_AUTOCOMPLETE]: SelectAutocompleteInput
  };

  let inputType =
    input.type && Object.values(TYPE).includes(input.type)
      ? input.type
      : 'default';

  let InputComponent = INPUTS[inputType];

  return (
    <div className='view_item'>
      <div className='item_name'>{input.label || input.name}</div>
      <div className='item_view'>
        <InputHelp input={input}>
          <InputComponent input={input} onChange={onChange}></InputComponent>
        </InputHelp>
      </div>
    </div>
  );
}
