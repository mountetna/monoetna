import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import BooleanInput from '../boolean';
import {InputSpecification} from '../input_types';

describe('BooleanInput', () => {
  let input: InputSpecification;
  let onChange: jest.Mock;

  function inputChecked(component: ReactWrapper): boolean {
    return (
      component.find('input').filterWhere((n) => n.prop('checked') || false)
        .length === 1
    );
  }

  beforeEach(() => {
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
      name: 'test-input',
      data: null
    };

    onChange = jest.fn();
  });

  it('returns boolean value when onChange called', () => {
    // Really just making sure we're capturing the right
    //   property out of the event.
    let component = mount(<BooleanInput input={input} onChange={onChange} />);

    expect(inputChecked(component)).toEqual(false);

    component
      .find('input')
      .simulate('change', {target: {checked: true}})
      .update();

    expect(onChange).toHaveBeenCalledWith('test-input', true);
  });

  it('sets checked state when input.value provided', () => {
    input.value = true;
    const component = mount(<BooleanInput input={input} onChange={onChange} />);

    expect(inputChecked(component)).toEqual(true);
  });
});
