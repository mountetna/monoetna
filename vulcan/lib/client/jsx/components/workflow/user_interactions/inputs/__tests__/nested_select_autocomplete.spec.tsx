import React from 'react';
import NestedSelectAutocompleteInput from '../nested_select_autocomplete';
import {DataEnvelope} from '../input_types';
import {integrateElement, setupBefore} from "../../../../../test_utils/integration";
import {Maybe, some} from "../../../../../selectors/maybe";
import {act, ReactTestInstance} from "react-test-renderer";
import {
  clickNode,
  includesClassNamePredicate,
  matchesTextPredicate, matchesTypePredicate,
  text
} from "../../../../../test_utils/rendered";

describe('NestedSelectAutocompleteInput', () => {
  const onChange = setupBefore(() => jest.fn());
  const value = setupBefore(() => null as Maybe<string>);
  const data = setupBefore(() => {
    return {
      'options-a': {
        option1: {
          suboption1: null,
          suboption2: {
            grandchild1: null,
            grandchild2: null
          }
        }
      },
      'options-b': {
        option2: {
          another1: {
            stepchild1: null,
            stepchild2: null
          }
        }
      }
    } as DataEnvelope<any>;
  });

  const integrated = setupBefore(() =>
    integrateElement(<NestedSelectAutocompleteInput
      onChange={onChange.value} value={value.value} data={data.value}/>))

  async function clickDropdownOption(component: ReactTestInstance, optionText: string) {
    await act(async () => {
      component
        .find(matchesTextPredicate(optionText))
        .props.onClick();
    })
  }

  async function clickIconWrapper(component: ReactTestInstance, idx: number) {
    await clickNode(component, includesClassNamePredicate('icon-wrapper'), idx);
  }

  async function clickLi(component: ReactTestInstance, idx: number) {
    await clickNode(component, matchesTypePredicate('li'), idx);
  }

  function getSelectedOption(component: ReactTestInstance) {
    return component.findAllByType('input')[0].props.value;
  }

  it('correctly manages child(ren) selects', async () => {
    const {node} = integrated.value;

    expect(node.root.findAllByType('input').length).toEqual(1);
    await clickIconWrapper(node.root, 0);
    await clickLi(node.root, 0);

    expect(getSelectedOption(node.root)).toEqual('option1');
    expect(node.root.findAllByType('input').length).toEqual(2);

    await clickIconWrapper(node.root, -1);

    expect(node.root.findAllByType('li').map((li) => text(li))).toEqual([
      'suboption1',
      'suboption2'
    ]);

    await clickIconWrapper(node.root, 0);
    // Close the second drop-down
    await clickIconWrapper(node.root, 1);
    // second item in the first drop-down
    await clickDropdownOption(node.root, 'option2');

    expect(node.root.findAllByType('input').length).toEqual(2);
    expect(getSelectedOption(node.root)).toEqual('option2');


    // Click the second drop-down again
    await clickIconWrapper(node.root, 1);

    expect(node.root.findAllByType('li').map(text)).toEqual(['another1']);
  });

  it('returns leaf value or null if not leaf', async () => {
    const {node} = integrated.value;
    expect(node.root.findAllByType('input').length).toEqual(1);
    await clickIconWrapper(node.root, 0);
    await clickLi(node.root, 0);
    expect(onChange.value).toBeCalledWith(null);

    await clickIconWrapper(node.root, -1);
    await clickDropdownOption(node.root, 'suboption1');
    expect(onChange.value).toBeCalledWith(some('suboption1'));

    await clickIconWrapper(node.root, -1);
    await clickDropdownOption(node.root, 'suboption2');

    expect(onChange.value).toBeCalledWith(null);
  });

  describe('with a given value', () => {
    value.replace(() => some('stepchild1'));

    it('can find an existing path when given a value', () => {
      const {node} = integrated.value;

      expect(node.root.findAllByType('input').length).toEqual(3);
      expect(node.root.findAllByType('input').map((i) => i.props.value)).toEqual([
        'option2',
        'another1',
        'stepchild1'
      ]);
    });
  })

  describe('with a stepchild set as a value', () => {
    value.replace(() => some('stepchild1'));

    it('correctly updates dropdowns when given a value', async () => {
      const {node} = integrated.value;

      expect(node.root.findAllByType('input').length).toEqual(3);
      await clickIconWrapper(node.root, 0);
      await clickLi(node.root, 0);
      expect(node.root.findAllByType('input').length).toEqual(2);
    });
  })
});
