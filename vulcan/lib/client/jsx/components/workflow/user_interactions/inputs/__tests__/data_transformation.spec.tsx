import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import DataTransformationInput from '../data_transformation';
import {BoundInputSpecification} from '../input_types';

describe('DataTransformationInput', () => {
  let input: BoundInputSpecification<string[], string[]>;
  let onChange: jest.Mock;

  beforeEach(() => {
    onChange = jest.fn();
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
      name: 'doesnotmatter',
      data: {
        col1: {0: 1, '1': 2, '2': 3},
        col2: {'0': 'abc', '1': '321', '2': 'xyz'}
      },
      onChange,
      source: ''
    };
  });

  it('mounts with all options checked', () => {
    const component = mount(
      <CheckboxesInput
        data={input.data}
        value={input.value}
        onChange={onChange}
      />
    );

    expect(component.find('.checkbox-input-option').length).toEqual(4);
    expect(renderedlabels(component)).toEqual(['1', '2', 'a', 'b']);
    expect(onChange).toHaveBeenCalledWith([['1', '2', 'a', 'b']]);
  });

  it('reflects input.value when provided', () => {
    input.value = [['2', 'b']];
    const component = mount(
      <CheckboxesInput
        data={input.data}
        value={input.value}
        onChange={onChange}
      />
    );

    expect(component.find('.checkbox-input-option').length).toEqual(4);
    expect(renderedlabels(component)).toEqual(['1', '2', 'a', 'b']);
    expect(onChange).not.toHaveBeenCalled();
    expect(checkedCheckboxesText(component)).toEqual(['2', 'b']);
  });
});
