import React from 'react';
import SingleDropdownMulticheckbox from '../single_dropdown_multicheckbox';
import {DataEnvelope} from '../input_types';
import {Maybe, some} from "../../../../../selectors/maybe";
import {integrateElement, setupBefore} from "../../../../../test_utils/integration";
import {
  findAllByClassName,
  includesClassNamePredicate,
  matchesTypePredicate, text
} from "../../../../../test_utils/rendered";
import {act, ReactTestInstance} from "react-test-renderer";

describe('SingleDropdownMulticheckbox', () => {
  const onChange = setupBefore(() => jest.fn());
  const value = setupBefore(() => null as Maybe<DataEnvelope<string[]>>);
  const data = setupBefore(() => {
    return {
      'options-a': {
        option1: ['1', '2', '3'],
        option2: ['x', 'y', 'z']
      },
      'options-b': {
        option3: ['9', '8', '7'],
        option4: ['a', 'b', 'c']
      }
    } as DataEnvelope<DataEnvelope<string[]>>
  });

  const integrated = setupBefore(() =>
    integrateElement(<SingleDropdownMulticheckbox
      onChange={onChange.value} value={value.value} data={data.value}/>))

  function clickCheckbox(component: ReactTestInstance, index: number) {
    component
      .findAll(includesClassNamePredicate('checkbox-input-option'))[index]
      .findByType('input')
      .props.onChange()
  }

  function renderedCheckboxesText(component: ReactTestInstance) {
    return component.findAll(
      includesClassNamePredicate('checkbox-input-option')).map(text);
  }

  function checkedCheckboxesText(component: ReactTestInstance) {
    return component.findAll(includesClassNamePredicate('checkbox-input-option'))
      .filter(node => node.findByType('input').props.checked)
      .map(text);
  }

  function renderedDropdownValue(component: ReactTestInstance) {
    return component
      .find(includesClassNamePredicate('dropdown-autocomplete-input'))
      .find(matchesTypePredicate('input'))
      .props.value;
  }
  it('correctly pre-selects all checkboxes if given `null` as value', async () => {
    const {node} = integrated.value;
    expect(findAllByClassName(node.root, 'dropdown-autocomplete-input').length).toEqual(1);
    expect(findAllByClassName(node.root, 'checkbox-input-option').length).toEqual(3);
    expect(renderedDropdownValue(node.root)).toEqual('option1');
    expect(renderedCheckboxesText(node.root)).toEqual(['1', '2', '3']);
    expect(onChange.value).toHaveBeenCalledWith(some({
      option1: ['1', '2', '3'],
      option2: ['x', 'y', 'z'],
      option3: ['9', '8', '7'],
      option4: ['a', 'b', 'c']
    }));
  });

  describe('with some selected options', () => {
    value.replace(() => {
      return some({
        option1: ['1', '2', '3'],
        option2: ['x', 'y', 'z'],
        option3: ['9', '8', '7'],
        option4: ['a', 'b', 'c']
      });
    })

    it('updates state when clicking checkboxes', () => {
      const {node} = integrated.value;

      clickCheckbox(node.root, 0);
      expect(onChange.value).toHaveBeenLastCalledWith(some({
        option1: ['2', '3'],
        option2: ['x', 'y', 'z'],
        option3: ['9', '8', '7'],
        option4: ['a', 'b', 'c']
      }));
    });

    describe('when missing some options in the value', () => {
      value.replace(() => {
        return some({
          option1: ['2'],
          option2: ['x', 'y', 'z'],
          option3: ['8', '7'],
          option4: ['a', 'b', 'c']
        });
      })

      it('updates state when clicking checkboxes', () => {
        const {node} = integrated.value;

        clickCheckbox(node.root, 0);
        expect(onChange.value).toHaveBeenLastCalledWith(some({
          option1: ['2', '1'],
          option2: ['x', 'y', 'z'],
          option3: ['8', '7'],
          option4: ['a', 'b', 'c']
        }));
      });
    })
  })

  it('checkboxes switch when dropdown value changes', async () => {
    const {node} = integrated.value;

    await act(async () => {
      node.root.find(includesClassNamePredicate('icon-wrapper'))
        .props.onClick()
    })

    const options = node.root.find(includesClassNamePredicate('dropdown-autocomplete-options'))
      .findAllByType('li')

    await act(async () => {
      options[options.length - 1].props.onClick();
    });

    expect(renderedDropdownValue(node.root)).toEqual('option4');
    expect(renderedCheckboxesText(node.root)).toEqual([
      'a',
      'b',
      'c'
    ]);
    expect(checkedCheckboxesText(node.root)).toEqual(['c']);
  });
});
