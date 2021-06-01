import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import CheckboxesInput from '../checkboxes';
import {InputSpecification} from '../input_types';

describe('CheckboxesInput', () => {
  let input: InputSpecification;
  let onChange: jest.Mock;

  function checkedCheckboxesText(component: ReactWrapper) {
    return component
      .find('.checkbox-input-option')
      .filterWhere((n) => n.find('input').prop('checked') || false)
      .map((n) => n.text());
  }

  function renderedlabels(component: ReactWrapper) {
    return component.find('.checkbox-input-option').map((n) => n.text());
  }

  beforeEach(() => {
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
      name: 'test-input',
      data: ['1', '2', 'a', 'b']
    };

    onChange = jest.fn();
  });

  it('mounts with all options checked', () => {
    const component = mount(
      <CheckboxesInput input={input} onChange={onChange} />
    );

    expect(component.find('.checkbox-input-option').length).toEqual(4);
    expect(renderedlabels(component)).toEqual(['1', '2', 'a', 'b']);

    expect(checkedCheckboxesText(component)).toEqual(['1', '2', 'a', 'b']);
  });
});
