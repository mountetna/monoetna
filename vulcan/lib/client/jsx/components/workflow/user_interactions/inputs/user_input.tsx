import React, {useContext, useEffect} from 'react';

import {defaultContext, VulcanContext} from '../../../../contexts/vulcan_context';
import {
  addValidationErrors,
  removeValidationErrors
} from '../../../../actions/vulcan_actions';
import InputHelp from './input_help';
import {
  InputBackendComponent,
  BoundInputSpecification,
  InputValidator
} from './input_types';
import { components, InputType, INPUT_TYPES } from '../../../ui_components';

function backendComponentOf(
  type: InputType
): [InputBackendComponent, InputValidator] {
  const result = type in components ? components[type] : components[INPUT_TYPES.STRING];
  return [result[0], result[1]];
}

export default function UserInput({
  input,
  hideLabel
}: {
  input: BoundInputSpecification;
  hideLabel?: boolean;
}) {
  const [InputComponent, Validator] = backendComponentOf(input.ui_component);
  const {dispatch} = useContext(VulcanContext);
  const {onChange, data, value, name, label, doc} = input;

  useEffect(() => {
    const errors = Validator({ data, value });
    if (errors.length > 0) {
      dispatch(
        addValidationErrors(name, label, errors)
      );
    }
    // On unmount, remove all errors associated with the input
    return () => {
      dispatch(removeValidationErrors(errors));
    };
  }, [dispatch, Validator, value, name, label, data]);

  return (
    <div className='view_item'>
      {!hideLabel ? (
        <div className='item_name'>{label}</div>
      ) : null}
      <div className='item_view'>
        <InputHelp doc={doc || ''}>
          {/* Input components should absolutely never have access to the top level context,
              as inputs are frequently nested in unpredictable ways and have no guarantee
              of their own association with the top level context.  Better to make input
              behavior stable by locking out the possibility of access to the greater context
              and ensuring the InputComponent interface is sufficiently detailed to capture
              its interface correctly in all cases.
           */}
          <VulcanContext.Provider value={defaultContext}>
            <InputComponent key={label}
              onChange={onChange}
              data={data}
              value={value}/>
          </VulcanContext.Provider>
        </InputHelp>
      </div>
    </div>
  );
}
