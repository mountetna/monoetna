import React, { useState } from 'react';
import NestedSelectAutocompleteInput from '../nested_select_autocomplete';
import {DataEnvelope} from '../input_types';
import {
  integrateElement,
  setupBefore
} from '../../../../../test_utils/integration';
import {Maybe, some} from '../../../../../selectors/maybe';
import {act, ReactTestInstance} from 'react-test-renderer';
import {matchesTypePredicate} from '../../../../../test_utils/rendered';
import Autocomplete from '@material-ui/lab/Autocomplete';

describe('NestedSelectAutocompleteInput', () => {
  const onChange = setupBefore(() => jest.fn());
  const defaultValue = setupBefore(() => null as Maybe<string | null>);
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
    integrateElement(() => {
      const [value, setValue] = useState(() => defaultValue.value);
      return <NestedSelectAutocompleteInput
        onChange={(v) => {
          setValue(v);
          onChange.value(v);
        }}
        value={value}
        data={data.value}
      />;
    })
  );

  function getSelectedOption(component: ReactTestInstance, index: number) {
    return component.findAll(matchesTypePredicate(Autocomplete))[index].props
      .value;
  }

  function getDropdownOptions(component: ReactTestInstance, index: number) {
    return component.findAll(matchesTypePredicate(Autocomplete))[index].props
      .options;
  }

  async function clickAutocompleteOption(
    component: ReactTestInstance,
    index: number,
    optionIndex: number
  ) {
    const autocomplete = component.findAll(matchesTypePredicate(Autocomplete))[
      index
    ];

    const options = autocomplete.props.options;
    // console.log("picking " + options[optionIndex])
    await act(async () => {
      autocomplete.props.onChange(null, options[optionIndex]);
    });
  }

  it('correctly manages child(ren) selects', async () => {
    const {node} = integrated.value;

    expect(node.root.findAllByType('input').length).toEqual(1);

    await clickAutocompleteOption(node.root, 0, 0);

    expect(getSelectedOption(node.root, 0)).toEqual('option1');
    expect(node.root.findAllByType('input').length).toEqual(2);

    expect(getDropdownOptions(node.root, 1)).toEqual([
      'suboption1',
      'suboption2'
    ]);

    await clickAutocompleteOption(node.root, 0, 1);

    expect(node.root.findAllByType('input').length).toEqual(2);
    expect(getSelectedOption(node.root, 0)).toEqual('option2');

    expect(getDropdownOptions(node.root, 1)).toEqual(['another1']);
  });

  it('returns leaf value or null if not leaf', async () => {
    const {node} = integrated.value;
    expect(node.root.findAllByType('input').length).toEqual(1);
    await clickAutocompleteOption(node.root, 0, 0);
    expect(onChange.value).toBeCalledWith(null);
    expect(node.root.findAllByType('input').length).toEqual(2);

    await clickAutocompleteOption(node.root, 1, 0);
    expect(node.root.findAllByType('input').length).toEqual(2);
    expect(onChange.value).toBeCalledWith(some('suboption1'));

    await clickAutocompleteOption(node.root, 1, 1);
    expect(node.root.findAllByType('input').length).toEqual(3);
    expect(onChange.value).toBeCalledWith(null);
  });

  describe('with a given value', () => {
    defaultValue.replace(() => some('stepchild1'));

    it('can find an existing path when given a value', async () => {
      const {node} = integrated.value;

      expect(node.root.findAllByType('input').length).toEqual(3);
      expect(
        node.root.findAllByType('input').map((i) => i.props.value)
      ).toEqual(['option2', 'another1', 'stepchild1']);
    });
  });

  describe('with a stepchild set as a value', () => {
    defaultValue.replace(() => some('stepchild1'));

    it('correctly updates dropdowns when given a value', async () => {
      let {node} = integrated.value;

      expect(node.root.findAllByType('input').length).toEqual(3);
      await clickAutocompleteOption(node.root, 0, 0);
      expect(node.root.findAllByType('input').length).toEqual(2);
    });
  });
});

// Methods not actually needed, but two ways to wait until renders complete if the issue is ever that await isn't actually waiting long enough.
// function delay(ms: number) {
//   return new Promise(resolve => setTimeout(resolve, ms));
// }
// while (true) {
//   await act(() => {})
//   if (node.root.findByType(NestedSelectAutocompleteInput).props.value == [null]) {
//     break
//   }
// }