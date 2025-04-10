import React, {useContext, useEffect, PropsWithChildren} from 'react';

import {defaultContext, VulcanContext} from '../../../contexts/vulcan_context';
import {
  addValidationErrors,
  removeValidationErrors
} from '../../../actions/vulcan_actions';
import {
  InputBackendComponent,
  BoundInputSpecification,
  InputValidator
} from '../ui_definitions/input_types';
import { inputComponents, InputType } from '../../ui_components';
import Icon from 'etna-js/components/icon';

function backendComponentOf(
  type: InputType
): [InputBackendComponent, InputValidator] {
  const result = type in inputComponents ? inputComponents[type] : inputComponents['default'];
  return [result[0], result[1]];
}

export function InputHelp({doc = '', children}: PropsWithChildren<{doc?: string}>) {
  return (
    <div className='input-help'>
      <div className='input-help-children-wrapper'>{children}</div>
      <div className='help-icon-wrapper' title={doc}>
        {doc ? (
          <Icon icon='question-circle' className='help-icon'/>
        ) : null}
      </div>
    </div>
  );
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
  const {onChange, data, value, name, label, doc, valueKeyMap, defaultValue} = input;

  useEffect(() => {
    const firstValueKey = Object.keys(valueKeyMap)[0]
    const errors = Validator({ data, value: value[firstValueKey] });
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
              value={value}
              defaultValue={defaultValue}/>
          </VulcanContext.Provider>
        </InputHelp>
      </div>
    </div>
  );
}
