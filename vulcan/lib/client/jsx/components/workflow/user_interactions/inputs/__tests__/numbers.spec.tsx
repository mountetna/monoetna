import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import {BoundInputSpecification} from '../input_types';
import { BaseTextFieldProps, FormHelperText, TextField } from '@material-ui/core';
import FloatInput from '../float';
import IntegerInput from '../integer';

function displayedValue(component: ReactWrapper) {
  return component.find(TextField).prop('value') as BaseTextFieldProps
}

function hasError(component: ReactWrapper) {
  return component.find(TextField).prop('error')
}

function callOnChange(component: ReactWrapper, args: any) {
  function call(component: ReactWrapper) {return component.find(TextField).invoke('onChange') as Function}
  return call(component)(args)
}

describe('FloatInput', () => {
  let input: BoundInputSpecification<number|null, number|null>;
  let onChange: jest.Mock;

  beforeEach(() => {
    onChange = jest.fn();
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
      name: 'doesnotmatter',
      data: { null: null },
      onChange,
      source: '',
    };
  });

  it('mounts with 0 when not given a value, and a label (label given)', () => {
    const component = mount(
      <FloatInput data={input.data} value={null} onChange={onChange} label={input.label}/>
    );

    expect(component.find(TextField).length).toEqual(1);
    expect(onChange).toHaveBeenCalledWith([0]); // = how 0 gets there
    expect(displayedValue(component)).toEqual('0')
    expect(hasError(component)).toEqual(false)

    expect(component.text()).toEqual('Abcdef');
  });

  it('reflects input.value when provided, and no label (not given)', () => {    
    const component = mount(
      <FloatInput data={input.data} value={[1]} onChange={onChange} />
    );

    expect(component.find(TextField).length).toEqual(1);
    expect(displayedValue(component)).toEqual('1')
    expect(hasError(component)).toEqual(false)

    expect(component.text()).toEqual('');
  });

  it('allows entry of non-numbers, but enters error state and call onChange with null', () => {
    const component = mount(
      <FloatInput data={input.data} value={input.value} onChange={onChange} />
    );

    // At start
    expect(displayedValue(component)).toEqual('0')
    expect(onChange).toHaveBeenCalledWith([0]);
    expect(hasError(component)).toEqual(false)

    // Given '1.5'
    callOnChange(component,{target: {value: '1.5'}})
    expect(displayedValue(component)).toEqual('1.5')
    expect(onChange).toHaveBeenCalledWith([1.5]);
    expect(hasError(component)).toEqual(false)

    // Given '1a'
    callOnChange(component,{target: {value: '1a'}})
    expect(displayedValue(component)).toEqual('1a')
    expect(onChange).toHaveBeenCalledWith([null]);
    expect(hasError(component)).toEqual(true)
  })
});

describe('IntegerInput', () => {
  let input: BoundInputSpecification<number|null, number|null>;
  let onChange: jest.Mock;

  beforeEach(() => {
    onChange = jest.fn();
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
      name: 'doesnotmatter',
      data: { null: null  },
      onChange,
      source: '',
    };
  });

  it('mounts with 0 when not given a value, and a label (label given)', () => {
    const component = mount(
      <IntegerInput data={input.data} value={null} onChange={onChange} label={input.label}/>
    );

    expect(component.find(TextField).length).toEqual(1);
    expect(onChange).toHaveBeenCalledWith([0]); // = how 0 gets there
    expect(displayedValue(component)).toEqual('0')
    expect(hasError(component)).toEqual(false)

    expect(component.text()).toEqual('Abcdef');
  });

  it('reflects input.value when provided, and no label (not given)', () => {    
    const component = mount(
      <IntegerInput data={input.data} value={[1]} onChange={onChange} />
    );

    expect(component.find(TextField).length).toEqual(1);
    expect(displayedValue(component)).toEqual('1')
    expect(hasError(component)).toEqual(false)

    expect(component.text()).toEqual('');
  });

  it('allows entry of non-numbers & non-integers, but enters error state and call onChange with null', () => {
    const component = mount(
      <IntegerInput data={input.data} value={input.value} onChange={onChange} />
    );

    // At start
    expect(displayedValue(component)).toEqual('0')
    expect(onChange).toHaveBeenCalledWith([0]);
    expect(hasError(component)).toEqual(false)

    // Given '1.5'
    callOnChange(component,{target: {value: '1.5'}})
    expect(displayedValue(component)).toEqual('1.5')
    expect(onChange).toHaveBeenCalledWith([null]);
    expect(hasError(component)).toEqual(true)

    // Given '1a'
    callOnChange(component,{target: {value: '1a'}})
    expect(displayedValue(component)).toEqual('1a')
    expect(onChange).toHaveBeenCalledWith([null]);
    expect(hasError(component)).toEqual(true)
  })
});
