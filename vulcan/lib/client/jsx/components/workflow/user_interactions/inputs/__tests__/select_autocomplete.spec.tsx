import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import SelectAutocompleteInput from '../select_autocomplete';
import {BoundInputSpecification} from '../input_types';
import { BaseTextFieldProps, ButtonBase, Checkbox, FormControlLabel, FormHelperText, Popper, TextField } from '@material-ui/core';
import Autocomplete from '@material-ui/lab/Autocomplete';
import { StringOptions } from '../monoids';

describe('SelectAutocompleteInput', () => {
  let input: BoundInputSpecification<string|null, StringOptions>;
  let onChange: jest.Mock;

  // Check options.filtered is shorter than options
  // Check options.display is length 1 if set maxOptions to 1
  // Check that helperText is generated when display < maxOptions
  // Check that value is null when nothing selected, but proper once an option is selected
  // Check that onChange can be overwritten

  beforeEach(() => {
    onChange = jest.fn();
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
      name: 'doesnotmatter',
      data: { 'options': ['1', '2', 'a', 'b'] },
      onChange,
      source: '',
    };
  });

  it('mounts with value null and a label (label given)', () => {
    const component = mount(
      <SelectAutocompleteInput data={input.data} value={input.value} onChange={onChange} label={input.label}/>
    );

    expect(component.find(Autocomplete).length).toEqual(1);
    expect(onChange).toHaveBeenCalledWith([null]);


    expect(component.text()).toEqual('Abcdef');
  });

  it('reflects input.value when provided, and no label (not given)', () => {    
    input.value = ['1']
    const component = mount(
      <SelectAutocompleteInput data={input.data} value={input.value} onChange={onChange} />
    );

    const textFieldProps = component.find(TextField).prop('inputProps') as BaseTextFieldProps
    expect(component.find(Autocomplete).length).toEqual(1);
    expect(onChange).not.toHaveBeenCalled();
    expect(textFieldProps['value']).toEqual('1')

    expect(component.text()).toEqual('');
  });

  it('truncates options based on maxOptions', () => {
    const component = mount(
      <SelectAutocompleteInput data={input.data} value={input.value} onChange={onChange} maxOptions={1}/>
    );
    expect(component.find(FormHelperText).length).toEqual(1)
    expect(component.text()).toEqual('4 options, only 1 shown')
  })
});
