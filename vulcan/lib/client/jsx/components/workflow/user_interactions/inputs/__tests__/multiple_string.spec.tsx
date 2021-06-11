import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import MultipleStringInput from '../multiple_string';
import {InputSpecification} from '../input_types';

describe('MultipleStringInput', () => {
  let input: InputSpecification;
  let onChange: jest.Mock;

  function renderedLabels(component: ReactWrapper) {
    return component.find('TextInput').map((n) => n.prop('header'));
  }

  function renderedValues(component: ReactWrapper) {
    return component.find('TextInput').map((n) => n.prop('value'));
  }

  beforeEach(() => {
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
      name: 'test-input',
      data: {
        'options-a': {
          1: '1',
          2: '2'
        },
        'options-b': {
          3: '3',
          4: '4'
        }
      }
    };

    onChange = jest.fn();
  });

  it('renders, initially', () => {
    const component = mount(
      <MultipleStringInput input={input} onChange={onChange} />
    );

    expect(component.find('TextInput').length).toEqual(0);
  });

  /*
  ToDo: test more properly / natively that 'data' gets transformed and given to 'value' upon an initial useEffect hook.
  'enzyme' does not handle this currently
  */
  it('renders, with values', () => {
    input.value = {
      1: '1',
      2: '2',
      3: '3',
      4: '4'
    };
    const component = mount(
      <MultipleStringInput input={input} onChange={onChange} />
    );

    expect(component.find('TextInput').length).toEqual(4);

    expect(renderedLabels(component)).toEqual([
      '1', '2', '3', '4'
    ]);
    expect(renderedValues(component)).toEqual([
      '1', '2', '3', '4'
    ]);
  });

  it('correctly sets the values when text is edited', () => {
    input.value = {
      1: '1',
      2: '2',
      3: '3',
      4: '4'
    };
    let component = mount(
      <MultipleStringInput input={input} onChange={onChange} />
    );

    component
      .find('.ti-input')
      .first()
      .simulate('change', { target : { value: 'T-thing' } } )
    expect(onChange).toHaveBeenCalledWith('test-input', {
      1: 'T-thing',
      2: '2',
      3: '3',
      4: '4'
    });
  });

});
