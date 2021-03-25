import React from 'react';

import {TYPE} from '../../../../api_types';
import InputHelp from './input_help';
import BooleanInput from './boolean';
import FloatInput from './float';
import IntegerInput from './integer';
import MultiselectStringInput from './multiselect_string';
import SelectAutocompleteInput from './select_autocomplete';
import StringInput from './string';
import CheckboxesInput from './checkboxes';
import {InputBackendComponent, InputOnChange, InputSpecification, InputType} from "./types";

function backendComponentOf(type: InputType): InputBackendComponent {
  switch(type) {
    case TYPE.FLOAT:
      return FloatInput;
    case TYPE.INTEGER:
      return IntegerInput;
    case TYPE.BOOL:
      return BooleanInput;
    case TYPE.MULTISELECT_STRING:
      return MultiselectStringInput;
    case TYPE.SELECT_AUTOCOMPLETE:
      return SelectAutocompleteInput;
    case TYPE.CHECKBOXES:
      return CheckboxesInput;
    default:
      return StringInput;
  }
}

export default function UserInput({input, onChange, hideLabel}: {input: InputSpecification, onChange: InputOnChange, hideLabel: boolean}) {
  const InputComponent = backendComponentOf(input.type);

  return (
    <div className='view_item'>
      {!hideLabel ? (
        <div className='item_name'>{input.label || input.name}</div>
      ) : null}
      <div className='item_view'>
        <InputHelp input={input}>
          <InputComponent input={input} onChange={onChange}/>
        </InputHelp>
      </div>
    </div>
  );
}
