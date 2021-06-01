import React, {useContext, useEffect} from 'react';

import {VulcanContext} from '../../../../contexts/vulcan_context';
import {
  addValidationErrors,
  removeValidationErrors
} from '../../../../actions/vulcan';
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
  InputType,
  InputValidator
} from './input_types';
import NestedSelectAutocompleteInput from './nested_select_autocomplete';
import MultiselectStringAllInput from './multiselect_string_all';
import MultipleMultiselectStringAllInput from './multiple_multiselect_string_all';
import SimpleNullValidator from './validators/simple_null_validator';

function backendComponentOf(
  type: InputType
): [InputBackendComponent, InputValidator] {
  switch (type) {
    case TYPE.FLOAT:
      return [FloatInput, SimpleNullValidator];
    case TYPE.INTEGER:
      return [IntegerInput, SimpleNullValidator];
    case TYPE.BOOL:
      return [BooleanInput, SimpleNullValidator];
    case TYPE.MULTISELECT_STRING:
      return [MultiselectStringInput, SimpleNullValidator];
    case TYPE.SELECT_AUTOCOMPLETE:
      return [SelectAutocompleteInput, SimpleNullValidator];
    case TYPE.CHECKBOXES:
      return [CheckboxesInput, SimpleNullValidator];
    case TYPE.NESTED_SELECT_AUTOCOMPLETE:
      return [NestedSelectAutocompleteInput, SimpleNullValidator];
    case TYPE.MULTISELECT_STRING_ALL:
      return [MultiselectStringAllInput, SimpleNullValidator];
    case TYPE.MULTIPLE_MULTISELECT_STRING_ALL:
      return [MultipleMultiselectStringAllInput, SimpleNullValidator];
    default:
      return [StringInput, SimpleNullValidator];
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
  const [InputComponent, Validator] = backendComponentOf(input.type);
  const {dispatch} = useContext(VulcanContext);

  useEffect(() => {
    let errors = Validator(input);
    if (errors.length > 0) {
      dispatch(addValidationErrors(input.name, errors));
    } else {
      dispatch(removeValidationErrors(input.name));
    }
  }, [input.value]);

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
