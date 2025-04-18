import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import CheckboxesInput from '../checkboxes';
import {BoundInputSpecification} from '../../../input_types';
import { Checkbox, FormControlLabel } from '@material-ui/core';

describe('CheckboxesInput', () => {
  let input: BoundInputSpecification<string[], string[]>;
  let onChange: jest.Mock;
  let showError: jest.Mock;

  function checkedCheckboxesText(component: ReactWrapper) {
    return component
      .find(FormControlLabel)
      .filterWhere((n) => n.find('input').prop('checked') || false)
      .map((n) => n.text());
  }

  function renderedlabels(component: ReactWrapper) {
    return component.find(FormControlLabel).map((n) => n.text());
  }

  beforeEach(() => {
    onChange = jest.fn();
    showError = jest.fn();
    input = {
      value: {picked: null},
      label: 'Abcdef',
      name: 'doesnotmatter',
      data: { 'options': ['1', '2', 'a', 'b'] },
      defaultValue: ['1', '2', 'a', 'b'],
      onChange,
      showError,
      ui_component: 'doesnotmatter',
      valueKeyMap: {picked: 'doesnotmatter'}
    };
  });

  it('mounts with all options checked', () => {
    const component = mount(
      <CheckboxesInput data={input.data} value={input.value} onChange={onChange} />
    );

    expect(component.find(Checkbox).length).toEqual(4);
    expect(renderedlabels(component)).toEqual(['1', '2', 'a', 'b']);
    expect(onChange).toHaveBeenCalledWith({picked: [['1', '2', 'a', 'b']]});
  });

  it('reflects input.value when provided', () => {
    input.value = {picked: [['2', 'b']]};
    const component = mount(
      <CheckboxesInput data={input.data} value={input.value} onChange={onChange} />
    );

    expect(component.find(Checkbox).length).toEqual(4);
    expect(renderedlabels(component)).toEqual(['1', '2', 'a', 'b']);
    expect(onChange).not.toHaveBeenCalled();
    expect(checkedCheckboxesText(component)).toEqual(['2', 'b']);
  });
});
