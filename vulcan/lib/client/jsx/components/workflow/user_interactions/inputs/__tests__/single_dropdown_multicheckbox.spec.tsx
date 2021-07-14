import React from 'react';
import {Provider} from 'react-redux';
import {mount, ReactWrapper, ShallowWrapper} from 'enzyme';
import {mockStore} from 'etna-js/spec/helpers';
import SingleDropdownMulticheckbox from '../single_dropdown_multicheckbox';
import {InputSpecification} from '../input_types';
import {
  VulcanState,
  defaultVulcanState
} from '../../../../../reducers/vulcan_reducer';
import {VulcanProvider} from '../../../../../contexts/vulcan_context';

describe('SingleDropdownMulticheckbox', () => {
  let input: InputSpecification;
  let state: VulcanState;
  let onChange: jest.Mock;
  let store: any;

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

  function renderedDropdownValue(component: ReactWrapper) {
    return component
      .find('.dropdown-autocomplete-input')
      .first()
      .find('input')
      .prop('value');
  }

  beforeEach(() => {
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
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

    state = {...defaultVulcanState};

    store = mockStore({});
  });

  it('correctly pre-selects all checkboxes if given `null` as value', async () => {
    const component = mount(
      <SingleDropdownMulticheckbox input={input} onChange={onChange} />
    );

    expect(component.find('.dropdown-autocomplete-input').length).toEqual(1);
    expect(component.find('.checkbox-input-option').length).toEqual(3);
    expect(renderedDropdownValue(component)).toEqual('option1');
    expect(renderedCheckboxesText(component)).toEqual(['1', '2', '3']);
    expect(onChange).toHaveBeenCalledWith('test-input', {
      option1: ['1', '2', '3'],
      option2: ['x', 'y', 'z'],
      option3: ['9', '8', '7'],
      option4: ['a', 'b', 'c']
    });
  });

  it('updates state when clicking checkboxes', () => {
    input.value = {
      option1: ['1', '2', '3'],
      option2: ['x', 'y', 'z'],
      option3: ['9', '8', '7'],
      option4: ['a', 'b', 'c']
    };
    const component = mount(
      <SingleDropdownMulticheckbox input={input} onChange={onChange} />
    );

    expect(component.find('.dropdown-autocomplete-input').length).toEqual(1);
    expect(component.find('.checkbox-input-option').length).toEqual(3);
    expect(renderedDropdownValue(component)).toEqual('option1');
    expect(renderedCheckboxesText(component)).toEqual(['1', '2', '3']);
    expect(checkedCheckboxesText(component)).toEqual(['1', '2', '3']);

    clickCheckbox(component, 0);
    expect(onChange).toHaveBeenLastCalledWith('test-input', {
      option1: ['2', '3'],
      option2: ['x', 'y', 'z'],
      option3: ['9', '8', '7'],
      option4: ['a', 'b', 'c']
    });
  });

  it('correctly sets the checkbox states when given a default', () => {
    input.value = {
      option1: ['1', '2'],
      option2: ['y'],
      option3: ['8', '9', '7'],
      option4: ['c']
    };
    const component = mount(
      <SingleDropdownMulticheckbox input={input} onChange={onChange} />
    );

    expect(component.find('.dropdown-autocomplete-input').length).toEqual(1);
    expect(component.find('.checkbox-input-option').length).toEqual(3);
    expect(renderedDropdownValue(component)).toEqual('option1');
    expect(renderedCheckboxesText(component)).toEqual(['1', '2', '3']);
    expect(checkedCheckboxesText(component)).toEqual(['1', '2']);
  });

  it('checkboxes switch when dropdown value changes', async () => {
    input.value = {
      option1: ['1', '2'],
      option2: ['y'],
      option3: ['8', '9', '7'],
      option4: ['c']
    };

    const component = mount(
      <Provider store={store}>
        <VulcanProvider state={state}>
          <SingleDropdownMulticheckbox input={input} onChange={onChange} />
        </VulcanProvider>
      </Provider>
    );
    expect(component.find('.dropdown-autocomplete-input').length).toEqual(1);
    expect(component.find('.checkbox-input-option').length).toEqual(3);
    expect(renderedDropdownValue(component)).toEqual('option1');
    expect(renderedCheckboxesText(component)).toEqual(['1', '2', '3']);
    expect(checkedCheckboxesText(component)).toEqual(['1', '2']);

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

    let dropdownInput = '__test-input__dropdownValue';
    expect(renderedDropdownValue(component)).toEqual('option4');
    expect(onChange).toHaveBeenCalledWith(dropdownInput, 'option4');

    state = { ...state, bufferedInputValue: 'option4' };
    // Super hokey, but can't setState on functional components,
    //   to verify behavior. So we'll remount with the new state
    //   and check that the render is different
    const rerenderedComponent = mount(
      <Provider store={store}>
        <VulcanProvider state={state}>
          <SingleDropdownMulticheckbox input={input} onChange={onChange} />
        </VulcanProvider>
      </Provider>
    );

    expect(renderedCheckboxesText(rerenderedComponent)).toEqual([
      'a',
      'b',
      'c'
    ]);
    expect(checkedCheckboxesText(rerenderedComponent)).toEqual(['c']);
  });
});
