import React from 'react';
import {mount, ReactWrapper, ShallowWrapper} from 'enzyme';
import SingleDropdownMulticheckbox from '../single_dropdown_multicheckbox';
import {InputSpecification} from '../input_types';

describe('SingleDropdownMulticheckbox', () => {
  let input: InputSpecification;
  let onChange: jest.Mock;

  function clickCheckbox(component: ReactWrapper, index: number) {
    component
      .find('.checkbox-input-option')
      .at(index)
      .find('input')
      .simulate('change')
      .update();
  }

  function renderedCheckboxesText(component: ReactWrapper) {
    return component.find('.checkbox-input-option').map((n) => n.text());
  }

  function checkedCheckboxesText(component: ReactWrapper) {
    return component
      .find('.checkbox-input-option')
      .filterWhere((n) => n.find('input').prop('checked') || false)
      .map((n) => n.text());
  }

  function renderedDropdownValues(component: ReactWrapper) {
    return component
      .find('.dropdown-autocomplete-input')
      .map((n) => n.find('input').first().props().value);
  }

  beforeEach(() => {
    input = {
      type: 'doesnotmatter',
      default: null,
      label: 'Abcdef',
      name: 'test-input',
      data: {
        'options-a': {
          option1: ['1', '2', '3'],
          option2: ['x', 'y', 'z']
        },
        'options-b': {
          option3: ['9', '8', '7'],
          option4: ['a', 'b', 'c']
        }
      }
    };

    onChange = jest.fn();
  });

  it('correctly renders checkboxes as pre-selected', () => {
    const component = mount(
      <SingleDropdownMulticheckbox input={input} onChange={onChange} />
    );
    console.log(component);
    expect(component.find('.dropdown-autocomplete-input').length).toEqual(2);
    expect(component.find('.checkbox-input-option').length).toEqual(6);
    expect(renderedDropdownValues(component)).toEqual(['option1', 'option3']);
    expect(renderedCheckboxesText(component)).toEqual([
      '1',
      '2',
      '3',
      '7',
      '8',
      '9'
    ]);
  });

  it('updates state when clicking checkboxes', () => {
    const component = mount(
      <SingleDropdownMulticheckbox input={input} onChange={onChange} />
    );

    expect(component.find('.dropdown-autocomplete-input').length).toEqual(2);
    expect(component.find('.checkbox-input-option').length).toEqual(6);
    expect(renderedDropdownValues(component)).toEqual(['option1', 'option3']);
    expect(renderedCheckboxesText(component)).toEqual([
      '1',
      '2',
      '3',
      '7',
      '8',
      '9'
    ]);

    clickCheckbox(component, 0);
    expect(onChange).toHaveBeenLastCalledWith('test-input', {
      'options-a': {
        option1: ['2', '3'],
        option2: ['x', 'y', 'z']
      },
      'options-b': {
        option3: ['9', '8', '7'],
        option4: ['a', 'b', 'c']
      }
    });
  });

  it('resets the input in state when any set of checkboxes is []', () => {
    const component = mount(
      <SingleDropdownMulticheckbox input={input} onChange={onChange} />
    );
    expect(component.find('.dropdown-autocomplete-input').length).toEqual(2);
    expect(renderedDropdownValues(component)).toEqual(['option1', 'option3']);
    clickCheckbox(component, 0);
    expect(onChange).toHaveBeenCalledWith('test-input', {
      'options-a': {
        option1: ['2', '3'],
        option2: ['x', 'y', 'z']
      },
      'options-b': {
        option3: ['9', '8', '7'],
        option4: ['a', 'b', 'c']
      }
    });

    clickCheckbox(component, 1);
    expect(onChange).toHaveBeenCalledWith('test-input', {
      'options-a': {
        option1: ['3'],
        option2: ['x', 'y', 'z']
      },
      'options-b': {
        option3: ['9', '8', '7'],
        option4: ['a', 'b', 'c']
      }
    });
    clickCheckbox(component, 2);
    expect(onChange).toHaveBeenCalledWith('test-input', null);
  });

  it('correctly sets the checkbox states when given a default', () => {
    input.default = {
      'options-a': {
        option1: ['1', '2'],
        option2: ['y']
      },
      'options-b': {
        option3: ['8', '9', '7'],
        option4: ['c']
      }
    };
    const component = mount(
      <SingleDropdownMulticheckbox input={input} onChange={onChange} />
    );

    expect(component.find('.dropdown-autocomplete-input').length).toEqual(2);
    expect(component.find('.checkbox-input-option').length).toEqual(6);
    expect(renderedDropdownValues(component)).toEqual(['option1', 'option3']);
    expect(renderedCheckboxesText(component)).toEqual([
      '1',
      '2',
      '3',
      '7',
      '8',
      '9'
    ]);
    expect(checkedCheckboxesText(component)).toEqual(['1', '2', '7', '8', '9']);
  });

  it('checkboxes switch when dropdown value changes', () => {
    input.default = {
      'options-a': {
        option1: ['1', '2'],
        option2: ['y']
      },
      'options-b': {
        option3: ['8', '9', '7'],
        option4: ['c']
      }
    };
    const component = mount(
      <SingleDropdownMulticheckbox input={input} onChange={onChange} />
    );
    expect(component.find('.dropdown-autocomplete-input').length).toEqual(2);
    expect(component.find('.checkbox-input-option').length).toEqual(6);
    expect(renderedDropdownValues(component)).toEqual(['option1', 'option3']);
    expect(renderedCheckboxesText(component)).toEqual([
      '1',
      '2',
      '3',
      '7',
      '8',
      '9'
    ]);
    expect(checkedCheckboxesText(component)).toEqual(['1', '2', '7', '8', '9']);

    component
      .find('.icon-wrapper')
      .first()
      .simulate('click')
      .update()
      .find('.dropdown-autocomplete-options')
      .find('li')
      .last()
      .simulate('click')
      .update();

    expect(renderedDropdownValues(component)).toEqual(['option2', 'option3']);
    expect(renderedCheckboxesText(component)).toEqual([
      'x',
      'y',
      'z',
      '7',
      '8',
      '9'
    ]);
    expect(checkedCheckboxesText(component)).toEqual(['y', '7', '8', '9']);
  });
});
