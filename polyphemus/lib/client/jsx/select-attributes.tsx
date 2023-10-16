import React, { useContext } from 'react';
import DropdownAutocompleteMultiPickInput from 'etna-js/components/inputs/dropdown_autocomplete_mui_multi_choice';
import {MagmaContext} from 'etna-js/contexts/magma-context';

const SelectAttributes = ({value, update, modelName, filter, label, placeholder, reportError = false, disabled = false}:{
  value: string[];
  update: (e: string[]) => void;
  modelName: string;
  filter?: (a: any) => boolean;
  label?: string;
  placeholder?: string;
  reportError?: boolean;
  disabled?: boolean;
}) => {
  const {models} = useContext(MagmaContext);

  if (!label) label = `Select ${modelName} attributes`

  const attribute_names = Object.values(
    models[modelName]?.template?.attributes || {}
  )
    .filter((a: any) => !a.hidden && (!filter || filter(a)))
    .map((a: any) => a.name)
    .sort();

  return <DropdownAutocompleteMultiPickInput
    optionSet={attribute_names}
    value={value}
    label={label}
    placeholder={placeholder}
    disableClearable={false}
    disabled={disabled}
    enforceError={reportError}
    onChange={(event: any, e: string[]) => update(e)}
  />
};

export default SelectAttributes;
