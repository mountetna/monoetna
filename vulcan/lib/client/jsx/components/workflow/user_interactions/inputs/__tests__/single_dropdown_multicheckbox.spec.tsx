import React from 'react';

import SingleDropdownMulticheckbox from '../single_dropdown_multicheckbox';
import {DataEnvelope} from '../input_types';
import {Maybe, some} from '../../../../../selectors/maybe';
import {
  integrateElement,
  setupBefore
} from '../../../../../test_utils/integration';
import {
  includesClassNamePredicate,
  matchesTypePredicate,
  text
} from '../../../../../test_utils/rendered';
import {act, ReactTestInstance} from 'react-test-renderer';
import {Checkbox, FormControlLabel} from '@material-ui/core';
import Autocomplete from '@material-ui/lab/Autocomplete';

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
    } as DataEnvelope<DataEnvelope<string[]>>;
  });

  const integrated = setupBefore(() =>
    integrateElement(
      <SingleDropdownMulticheckbox
        onChange={onChange.value}
        value={value.value}
        data={data.value}
      />
    )
  );

  function clickCheckbox(component: ReactTestInstance, index: number) {
    // Autocomplete is index 0, so to get a checkbox at `index`, we'll add one.
    component.findAll(matchesTypePredicate('input'))[index + 1].props.onChange({
      target: {
        value: 'changed'
      }
    });
  }

  function renderedCheckboxesText(component: ReactTestInstance) {
    return component.findAll(matchesTypePredicate(FormControlLabel)).map(text);
  }

  function checkedCheckboxesText(component: ReactTestInstance) {
    return component
      .findAll(matchesTypePredicate(FormControlLabel))
      .filter((node) => node.findByType('input').props.checked)
      .map(text);
  }

  function renderedDropdownValue(component: ReactTestInstance) {
    return component
      .find(matchesTypePredicate(Autocomplete))
      .find(matchesTypePredicate('input')).props.value;
  }
  it('correctly pre-selects all checkboxes if given `null` as value', async () => {
    const {node} = integrated.value;
    expect(
      node.root.findAll(matchesTypePredicate(Autocomplete)).length
    ).toEqual(1);
    expect(node.root.findAll(matchesTypePredicate(Checkbox)).length).toEqual(3);
    expect(renderedDropdownValue(node.root)).toEqual('option1');
    expect(renderedCheckboxesText(node.root)).toEqual(['1', '2', '3']);
    expect(onChange.value).toHaveBeenCalledWith(
      some({
        option1: ['1', '2', '3'],
        option2: ['x', 'y', 'z'],
        option3: ['9', '8', '7'],
        option4: ['a', 'b', 'c']
      })
    );
  });

  describe('with some selected options', () => {
    value.replace(() => {
      return some({
        option1: ['1', '2', '3'],
        option2: ['x', 'y', 'z'],
        option3: ['9', '8', '7'],
        option4: ['a', 'b', 'c']
      });
    });

    it('updates state when clicking checkboxes', () => {
      const {node} = integrated.value;

      clickCheckbox(node.root, 0);
      expect(onChange.value).toHaveBeenLastCalledWith(
        some({
          option1: ['2', '3'],
          option2: ['x', 'y', 'z'],
          option3: ['9', '8', '7'],
          option4: ['a', 'b', 'c']
        })
      );
    });

    describe('when missing some options in the value', () => {
      value.replace(() => {
        return some({
          option1: ['2'],
          option2: ['x', 'y', 'z'],
          option3: ['8', '7'],
          option4: ['a', 'b']
        });
      });

      it('checkboxes switch when dropdown value changes', async () => {
        const {node} = integrated.value;

        const autocomplete = node.root.find(matchesTypePredicate(Autocomplete));

        const options = autocomplete.props.options;

        expect(renderedDropdownValue(node.root)).toEqual('option1');

        await act(async () => {
          autocomplete.props.onChange(null, options[options.length - 1]);
        });

        expect(renderedDropdownValue(node.root)).toEqual('option4');
        expect(renderedCheckboxesText(node.root)).toEqual(['a', 'b', 'c']);
        expect(checkedCheckboxesText(node.root)).toEqual(['a', 'b']);
      });

      it('updates state when clicking checkboxes', () => {
        const {node} = integrated.value;

        clickCheckbox(node.root, 0);
        expect(onChange.value).toHaveBeenLastCalledWith(
          some({
            option1: ['2', '1'],
            option2: ['x', 'y', 'z'],
            option3: ['8', '7'],
            option4: ['a', 'b']
          })
        );
      });
    });
  });
});
