import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import MultiselectStringInput from '../multiselect_string';
import {Maybe} from '../../../../../../selectors/maybe';
import { DataEnvelope } from '../../../input_types';

describe('MultiselectStringInput', () => {
  let data: {options: string[]};
  let value: Maybe<DataEnvelope<string[]>>;
  let onChange: jest.Mock;

  function openDropdown(component: ReactWrapper) {
    component
      .find('.add_item')
      .first()
      .simulate('click')
      .update()
      .find('.icon-wrapper')
      .first()
      .simulate('click')
      .update();
  }

  function clickOption(
    component: ReactWrapper,
    optionIndex: number
  ) {
    openDropdown(component);
    component.find('li').at(optionIndex).simulate('click').update();
  }

  function renderedItemsList(component: ReactWrapper) {
    return component.find('.delete_link').map((n) => n.text());
  }

  beforeEach(() => {
    value = {picked: null};
    data = {
        'options': ['1', '2', '3', 'x', 'y', 'z']
    };

    onChange = jest.fn();
  });

  it('Correctly selects first choice', () => {
    const component = mount(
      <MultiselectStringInput data={data} value={value} onChange={onChange} />
    );
    clickOption(component, 0);
    expect(onChange).toHaveBeenLastCalledWith({picked: [['1']]});
  });

  it('Correctly shows the list when given a value', () => {
    value = {
      picked: [['1', '2']]
    };

    const component = mount(
      <MultiselectStringInput data={data} value={value} onChange={onChange} />
    );

    expect(renderedItemsList(component)).toEqual([
      '1',
      '2'
    ]);
  });

  it('can remove a single entry with a value', () => {
    value = {
      picked: [['1', '2']]
    };

    const component = mount(
      <MultiselectStringInput data={data} value={value} onChange={onChange} />
    );

    component.find('.delete_link').last().simulate('click');
    expect(onChange).toHaveBeenLastCalledWith({picked: [['1']]});
  });
});
