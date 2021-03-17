import React from 'react';

import {TYPE} from '../../../../models/steps';
import InputHelp from './input_help';
import BooleanInput from './boolean';
import FloatInput from './float';
import IntegerInput from './integer';
import MultiselectStringInput from './multiselect_string';
import SelectAutocompleteInput from './select_autocomplete';
import StringInput from './string';
import {wrapPaneItem} from '../../../../utils/workflow';

export default function UserInput({input, onChange}) {
  const INPUTS = {
    default: StringInput,
    [TYPE.INTEGER]: IntegerInput,
    [TYPE.FLOAT]: FloatInput,
    [TYPE.BOOL]: BooleanInput,
    [TYPE.MULTISELECT_STRING]: MultiselectStringInput,
    [TYPE.SELECT_AUTOCOMPLETE]: SelectAutocompleteInput
  };

  let inputType = input.type || 'default';

  let Component = INPUTS[inputType];

  return wrapPaneItem(
    {
      name: input.name,
      value: (
        <InputHelp input={input}>
          <Component input={input} onChange={onChange}></Component>
        </InputHelp>
      )
    },
    input.name
  );
}
