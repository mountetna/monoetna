import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import MultipleMultiselectStringAllInput from '../multiple_multiselect_string_all';
import {InputSpecification} from '../input_types';

describe('MultipleMultiselectStringAllInput', () => {
  let input: InputSpecification;
  let onChange: jest.Mock;

  function addItem(
    component: ReactWrapper,
    inputIndex: number,
    optionIndex: number
  ) {
    component
      .find('.view_item')
      .at(inputIndex)
      .find('.add_item')
      .first()
      .simulate('click')
      .update()
      .find('.icon-wrapper')
      .first()
      .simulate('click')
      .update()
      .find('li')
      .at(optionIndex)
      .simulate('click')
      .update();
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

  it('correctly manages child(ren) selects', () => {
    const component = mount(
      <MultipleMultiselectStringAllInput input={input} onChange={onChange} />
    );
    expect(component.find('.view_item').length).toEqual(4);
    expect(component.find('.item_name').map((n) => n.text())).toEqual([
      'option1',
      'option2',
      'option3',
      'option4'
    ]);
    addItem(component, 0, 0);

    expect(component.find('.delete_link').first().text()).toEqual('1');
    expect(onChange).toHaveBeenCalledWith('test-input', null);

    addItem(component, 0, 1);
    expect(onChange).toHaveBeenCalledWith('test-input', null);

    addItem(component, 1, 2);
    expect(onChange).toHaveBeenCalledWith('test-input', null);

    addItem(component, 2, 1);
    expect(onChange).toHaveBeenCalledWith('test-input', null);
    addItem(component, 3, 2);

    expect(onChange).toHaveBeenCalledWith('test-input', null);
    addItem(component, 3, 0);

    expect(component.find('.delete_link').last().text()).toEqual('a');
    expect(onChange).toHaveBeenCalledWith('test-input', {
      'options-a': {
        option1: ['1', '2'],
        option2: ['z']
      },
      'options-b': {
        option3: ['8'],
        option4: ['c', 'a']
      }
    });
  });

  it('correctly sets the lists when given a default', () => {
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
      <MultipleMultiselectStringAllInput input={input} onChange={onChange} />
    );

    expect(component.find('.view_item').length).toEqual(4);
    expect(component.find('.item_name').map((n) => n.text())).toEqual([
      'option1',
      'option2',
      'option3',
      'option4'
    ]);
    expect(component.find('.delete_link').map((val) => val.text())).toEqual([
      '1',
      '2',
      'y',
      '8',
      '9',
      '7',
      'c'
    ]);
  });
});
