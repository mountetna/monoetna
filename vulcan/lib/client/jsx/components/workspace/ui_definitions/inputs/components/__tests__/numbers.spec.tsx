import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import {BoundInputSpecification} from '../../../input_types';
import TextField, { BaseTextFieldProps } from '@material-ui/core/TextField';
import FloatInput from '../float';
import IntegerInput from '../integer';

function displayedValue(component: ReactWrapper) {
  return component.find(TextField).prop('value') as BaseTextFieldProps;
}

function hasError(component: ReactWrapper) {
  return component.find(TextField).prop('error');
}

function callOnChange(component: ReactWrapper, args: any) {
  function call(component: ReactWrapper) {return component.find(TextField).invoke('onChange') as Function;}
  return call(component)(args);
}

describe('FloatInput', () => {
  let input: BoundInputSpecification<number|null, number|null>;
  let onChange: jest.Mock;
  let showError: jest.Mock;

  beforeEach(() => {
    onChange = jest.fn();
    showError = jest.fn();
    input = {
      value: null,
      label: 'Abcdef',
      name: 'doesnotmatter',
      data: { null: null },
      onChange,
      showError,
      ui_component: 'doesnotmatter',
      valueKeyMap: {picked: 'doesnotmatter'}
    };
  });

  it('mounts with 0.0 when not given a value, and a label (label given)', () => {
    const component = mount(
      <FloatInput data={input.data} value={{value: null}} onChange={onChange} label={input.label}/>
    );

    expect(component.find(TextField).length).toEqual(1);
    expect(onChange).toHaveBeenCalledWith({value: [0]}); // = how 0 gets there
    expect(displayedValue(component)).toEqual('0.0');

    expect(component.text()).toEqual('Abcdef');
  });

  it('reflects input.value when provided, and no label (not given)', () => {    
    const component = mount(
      <FloatInput data={input.data} value={{value: [1]}} onChange={onChange} />
    );

    expect(component.find(TextField).length).toEqual(1);
    expect(displayedValue(component)).toEqual('1.0');

    expect(component.text()).toEqual('');
  });

  it('does not allow entry of non-numbers', () => {
    const component = mount(
      <FloatInput data={input.data} value={{value: input.value}} onChange={onChange} />
    );

    // At start
    expect(displayedValue(component)).toEqual('0.0');
    expect(onChange).toHaveBeenCalledWith({value: [0]});

    // Given '1.5'
    callOnChange(component,{target: {value: '1.5'}});
    expect(onChange).toHaveBeenCalledWith({value: [1.5]});
    // ToFix: Not updating for the test?
    // expect(displayedValue(component)).toEqual('1.5');

    // Given '1a'
    callOnChange(component,{target: {value: '1a'}});
    expect(onChange).toHaveBeenLastCalledWith({value: [1.5]});
  });
});

describe('IntegerInput', () => {
  let input: BoundInputSpecification<number|null, number|null>;
  let onChange: jest.Mock;
  let showError: jest.Mock;

  beforeEach(() => {
    onChange = jest.fn();
    showError = jest.fn();
    input = {
      value: null,
      label: 'Abcdef',
      name: 'doesnotmatter',
      data: { data: null  },
      onChange,
      showError,
      ui_component: 'doesnotmatter',
      valueKeyMap: {picked: 'doesnotmatter'}
    };
  });

  it('mounts with 0 when not given a value, and a label (label given)', () => {
    const component = mount(
      <IntegerInput data={input.data} value={{value: null}} onChange={onChange} label={input.label}/>
    );

    expect(component.find(TextField).length).toEqual(1);
    expect(onChange).toHaveBeenCalledWith({value: [0]}); // = how 0 gets there
    expect(displayedValue(component)).toEqual('0');

    expect(component.text()).toEqual('Abcdef');
  });

  it('reflects input.value when provided, and no label (not given)', () => {
    const component = mount(
      <IntegerInput data={input.data} value={{value: [1]}} onChange={onChange} />
    );

    expect(component.find(TextField).length).toEqual(1);
    expect(displayedValue(component)).toEqual('1');

    expect(component.text()).toEqual('');
  });

  it('does not allow entry of non-numbers & non-integers', () => {
    const component = mount(
      <IntegerInput data={input.data} value={{value: input.value}} onChange={onChange} />
    );

    // At start
    expect(displayedValue(component)).toEqual('0');
    expect(onChange).toHaveBeenCalledWith({value: [0]});

    // Given '1.5'
    callOnChange(component,{target: {value: '1.5'}});
    expect(displayedValue(component)).toEqual('0');
    expect(onChange).toHaveBeenLastCalledWith({value: [0]});

    // Given '1a'
    callOnChange(component,{target: {value: '1a'}});
    expect(displayedValue(component)).toEqual('0');
    expect(onChange).toHaveBeenLastCalledWith({value: [0]});
  });
});
