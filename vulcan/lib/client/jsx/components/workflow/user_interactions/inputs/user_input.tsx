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
import {
  InputBackendComponent,
  InputOnChange,
  InputSpecification,
  InputType
} from './input_types';
import NestedSelectAutocompleteInput from './nested_select_autocomplete';
import MultiselectStringAllInput from './multiselect_string_all';
import MultipleMultiselectStringAllInput from './multiple_multiselect_string_all';

function backendComponentOf(type: InputType): InputBackendComponent {
  switch (type) {
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
    case TYPE.NESTED_SELECT_AUTOCOMPLETE:
      return NestedSelectAutocompleteInput;
    case TYPE.MULTISELECT_STRING_ALL:
      return MultiselectStringAllInput;
    case TYPE.MULTIPLE_MULTISELECT_STRING_ALL:
      return MultipleMultiselectStringAllInput;
    default:
      return StringInput;
  }
}

export default function UserInput({
  input,
  onChange,
  hideLabel
}: {
  input: InputSpecification;
  onChange: InputOnChange;
  hideLabel?: boolean;
}) {
  const InputComponent = backendComponentOf(input.type);

  return (
    <div className='view_item'>
      {!hideLabel ? (
        <div className='item_name'>{input.label || input.name}</div>
      ) : null}
      <div className='item_view'>
        <InputHelp input={input}>
          <InputComponent key={input.name} input={input} onChange={onChange} />
        </InputHelp>
      </div>
    </div>
  );
}
