import React, {useContext, useEffect} from 'react';

import {defaultContext, VulcanContext} from '../../../../contexts/vulcan_context';
import {
  addValidationErrors,
  removeValidationErrors
} from '../../../../actions/vulcan_actions';
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
import MultipleStringInput from './multiple_string';
import SingleDropdownMulticheckbox from './single_dropdown_multicheckbox';
import NotEmptyValidator from './validators/not_empty_validator';
import {
  AllInnerValuesNotEmptyValidator,
  AllInnerValuesNotEmptyValidatorStrong
} from './validators/all_inner_values_not_empty_validator';

function backendComponentOf(
  type: InputType
): [InputBackendComponent, InputValidator] {
  switch (type) {
    case TYPE.FLOAT:
      return [FloatInput, NotEmptyValidator];
    case TYPE.INTEGER:
      return [IntegerInput, NotEmptyValidator];
    case TYPE.BOOL:
      return [BooleanInput, NotEmptyValidator];
    case TYPE.MULTISELECT_STRING:
      return [MultiselectStringInput, NotEmptyValidator];
    case TYPE.SELECT_AUTOCOMPLETE:
      return [SelectAutocompleteInput, NotEmptyValidator];
    case TYPE.CHECKBOXES:
      return [CheckboxesInput, NotEmptyValidator];
    case TYPE.NESTED_SELECT_AUTOCOMPLETE:
      return [NestedSelectAutocompleteInput, NotEmptyValidator];
    case TYPE.MULTISELECT_STRING_ALL:
      return [MultiselectStringAllInput, NotEmptyValidator];
    case TYPE.MULTIPLE_MULTISELECT_STRING_ALL:
      return [
        MultipleMultiselectStringAllInput,
        AllInnerValuesNotEmptyValidator
      ];
    case TYPE.MULTIPLE_STRING:
      return [MultipleStringInput, AllInnerValuesNotEmptyValidatorStrong];
    case TYPE.SINGLE_DROPDOWN_MULTICHECKBOX:
      return [SingleDropdownMulticheckbox, AllInnerValuesNotEmptyValidator];
    default:
      return [StringInput, NotEmptyValidator];
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
      dispatch(
        addValidationErrors(input.name, input.label || input.name, errors)
      );
    } else {
      dispatch(removeValidationErrors(input.name));
    }
    // On unmount, remove all errors associated with the input
    return () => dispatch(removeValidationErrors(input.name));
  }, [input.value, input, dispatch, Validator]);

  return (
    <div className='view_item'>
      {!hideLabel ? (
        <div className='item_name'>{input.label || input.name}</div>
      ) : null}
      <div className='item_view'>
        <InputHelp input={input}>
          {/* Input components should absolutely never have access to the top level context,
              as inputs are frequently nested in unpredictable ways and have no guarantee
              of their own association with the top level context.  Better to make input
              behavior stable by locking out the possibility of access to the greater context
              and ensuring the InputComponent interface is sufficiently detailed to capture
              its interface correctly in all cases.
           */}
          <VulcanContext.Provider value={defaultContext}>
            <InputComponent key={input.name} input={input} onChange={onChange} />
          </VulcanContext.Provider>
        </InputHelp>
      </div>
    </div>
  );
}
