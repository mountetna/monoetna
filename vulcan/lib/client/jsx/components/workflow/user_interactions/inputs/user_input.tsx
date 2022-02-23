import React, {Dispatch, useContext, useEffect} from 'react';

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
  BoundInputSpecification,
  InputType,
  InputValidator
} from './input_types';
import NestedSelectAutocompleteInput from './nested_select_autocomplete';
import {ScatterPlotly, BarPlotly, YPlotly} from './visualizations';
import {
  NotEmptyValidator,
  StronglyNotEmptyValidator
} from './validators/not_empty_validator';
import {
  AllInnerValuesNotEmptyValidator,
} from './validators/all_inner_values_not_empty_validator';
import MultipleInput from "./multiple_input";
import SingleDropdownMulticheckbox from "./single_dropdown_multicheckbox";
import {stepOfSource} from "../../../../selectors/workflow_selectors";
import AllOutputValuesNotEmptyValidator from "./validators/all_output_values_not_empty_validator"
import DiffExpSC from './scDGE';
import DataTransformation from './data_transformation';
import AllInnerKeysNotNullValidator from './validators/all_inner_keys_not_null_validator';
import PlusSubsetValidator from './validators/PlusSubsetValidator';

const components: {[k: string]: [InputBackendComponent<any, any, any>, InputValidator<any, any>]} = {};
function configureComponent<Value, DataElement>(
  type: InputType,
  Input: InputBackendComponent<any, Value, DataElement>,
  validator: InputValidator<Value, DataElement>,
) {
  if (type in components) throw new Error(`Duplicate definition for ${type}`)
  components[type] = [Input, validator];
}

configureComponent(TYPE.STRING, StringInput, NotEmptyValidator);
configureComponent(TYPE.FLOAT, FloatInput, NotEmptyValidator);
configureComponent(TYPE.INTEGER, IntegerInput, NotEmptyValidator);
configureComponent(TYPE.BOOL, BooleanInput, NotEmptyValidator);
configureComponent(TYPE.SELECT_AUTOCOMPLETE, SelectAutocompleteInput, StronglyNotEmptyValidator);
configureComponent(TYPE.CHECKBOXES, CheckboxesInput, NotEmptyValidator);
configureComponent(TYPE.NESTED_SELECT_AUTOCOMPLETE, NestedSelectAutocompleteInput, StronglyNotEmptyValidator);
configureComponent(TYPE.MULTISELECT_STRING, MultiselectStringInput, NotEmptyValidator);
configureComponent(TYPE.MULTIPLE_STRING, MultipleInput(StringInput), AllInnerValuesNotEmptyValidator);
configureComponent(TYPE.DATA_TRANSFORMATION, DataTransformation, AllInnerKeysNotNullValidator);
configureComponent(TYPE.SCATTER_PLOTLY, ScatterPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator));
configureComponent(TYPE.BAR_PLOTLY, BarPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator));
configureComponent(TYPE.Y_PLOTLY, YPlotly, PlusSubsetValidator('rows_use',AllOutputValuesNotEmptyValidator));
configureComponent(TYPE.DIFF_EXP_SC, DiffExpSC, PlusSubsetValidator('subset',AllOutputValuesNotEmptyValidator));

configureComponent(TYPE.SINGLE_DROPDOWN_MULTICHECKBOX, SingleDropdownMulticheckbox, NotEmptyValidator);
configureComponent(TYPE.MULTIPLE_MULTISELECT_STRING_ALL, MultipleInput(MultiselectStringInput), AllInnerValuesNotEmptyValidator)


function backendComponentOf(
  type: InputType
): [InputBackendComponent, InputValidator] {
  const result = components[type];
  if (!result) return components[TYPE.STRING];
  return result;
}

export default function UserInput({
  input,
  hideLabel
}: {
  input: BoundInputSpecification;
  hideLabel?: boolean;
}) {
  const [InputComponent, Validator] = backendComponentOf(input.type);
  const {dispatch} = useContext(VulcanContext);
  const {onChange, data, value, source, label} = input;

  useEffect(() => {
    const errors = Validator({ data, value });
    if (errors.length > 0) {
      dispatch(
        addValidationErrors(stepOfSource(source) || null, label, errors)
      );
    }
    // On unmount, remove all errors associated with the input
    return () => {
      dispatch(removeValidationErrors(errors))
    };
  }, [dispatch, Validator, value, source, label, data]);
  console.log('input', input);
  return (
    <div className='view_item'>
      {!hideLabel ? (
        <div className='item_name'>{input.label}</div>
      ) : null}
      <div className='item_view'>
        <InputHelp doc={input.doc || ''}>
          {/* Input components should absolutely never have access to the top level context,
              as inputs are frequently nested in unpredictable ways and have no guarantee
              of their own association with the top level context.  Better to make input
              behavior stable by locking out the possibility of access to the greater context
              and ensuring the InputComponent interface is sufficiently detailed to capture
              its interface correctly in all cases.
           */}
          <VulcanContext.Provider value={defaultContext}>
            <InputComponent key={input.label} onChange={onChange} data={data} value={value} />
          </VulcanContext.Provider>
        </InputHelp>
      </div>
    </div>
  );
}
