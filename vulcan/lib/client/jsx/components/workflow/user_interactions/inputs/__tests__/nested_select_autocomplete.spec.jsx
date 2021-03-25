import React from 'react';
import renderer, {act} from 'react-test-renderer';
import {mount} from 'enzyme';
import NestedSelectAutocompleteInput from '../nested_select_autocomplete';

describe('NestedSelectAutocompleteInput', () => {
  let input;
  let onChange;

  function clickDropdownOption(component, optionText) {
    component
      .findWhere((n) => n.text() == optionText)
      .first()
      .simulate('click');
  }

  beforeEach(() => {
    input = {
      name: 'test-input',
      options: [
        {
          option1: {
            suboption1: null,
            suboption2: {
              grandchild1: null,
              grandchild2: null
            }
          },
          option2: {
            another1: {
              stepchild1: null,
              stepchild2: null
            }
          }
        }
      ]
    };

    onChange = jest.fn();
  });

  it('correctly manages child(ren) selects', () => {
    const component = mount(
      <NestedSelectAutocompleteInput
        input={input}
        onChange={onChange}
      ></NestedSelectAutocompleteInput>
    );

    expect(component.find('input').length).toEqual(1);

    component.find('.icon-wrapper').first().simulate('click');

    component.update();

    component.find('li').first().simulate('click');

    component.update();

    expect(component.find('input').first().instance().value).toEqual('option1');
    expect(component.find('input').length).toEqual(2);

    component.find('.icon-wrapper').last().simulate('click');

    component.update();

    expect(component.find('li').map((li) => li.text())).toEqual([
      'suboption1',
      'suboption2'
    ]);

    component.find('.icon-wrapper').first().simulate('click');

    // Close the second drop-down
    component.find('.icon-wrapper').at(1).simulate('click');

    component.update();

    // second item in the first drop-down
    clickDropdownOption(component, 'option2');

    component.update();

    expect(component.find('input').length).toEqual(2);
    expect(component.find('input').first().instance().value).toEqual('option2');

    component.find('.icon-wrapper').at(1).simulate('click');

    component.update();

    expect(component.find('li').map((li) => li.text())).toEqual(['another1']);
  });

  it('returns leaf value or null if not leaf', () => {
    const component = mount(
      <NestedSelectAutocompleteInput
        input={input}
        onChange={onChange}
      ></NestedSelectAutocompleteInput>
    );

    expect(component.find('input').length).toEqual(1);

    component.find('.icon-wrapper').first().simulate('click');

    component.update();

    component.find('li').first().simulate('click');

    component.update();

    expect(onChange).toBeCalledWith('test-input', null);

    component.find('.icon-wrapper').last().simulate('click');

    component.update();

    clickDropdownOption(component, 'suboption1');

    expect(onChange).toBeCalledWith('test-input', 'suboption1');

    component.find('.icon-wrapper').last().simulate('click');

    component.update();

    clickDropdownOption(component, 'suboption2');

    expect(onChange).toBeCalledWith('test-input', null);
  });

  fit('can find an existing path when given a default', () => {
    input.default = 'stepchild1';

    const component = mount(
      <NestedSelectAutocompleteInput
        input={input}
        onChange={onChange}
      ></NestedSelectAutocompleteInput>
    );

    expect(component.find('input').length).toEqual(3);
    expect(component.find('input').map((i) => i.instance().value)).toEqual([
      'option2',
      'another1',
      'stepchild1'
    ]);
  });
});
