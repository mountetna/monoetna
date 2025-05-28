import React, { useState } from 'react';
import {mount, ReactWrapper} from 'enzyme';
import DropdownInput from '../dropdown';
import {BoundInputSpecification} from '../../../input_types';
import Autocomplete from '@material-ui/lab/Autocomplete'
import TextField from '@material-ui/core/TextField';
import {act, ReactTestInstance} from 'react-test-renderer';
import {
  integrateElement
} from '../../../../../../test_utils/integration';
import {matchesTypePredicate} from '../../../../../../test_utils/rendered';

describe('DropdownInput', () => {
  let input: BoundInputSpecification<string[], string[]>;
  let onChange: jest.Mock;
  let showError: jest.Mock;
  let integrated: ReturnType<typeof integrateElement>;

  function initialValue(component: ReactWrapper) {
    return component
      .find(Autocomplete)
      .first()
      .prop('value');
  }

  function initiallabel(component: ReactWrapper) {
    return component.find(TextField).first().prop('label');
  }

  function getInitialDropdownOptions(component: ReactWrapper) {
    return component.find(Autocomplete).prop('options');
  }

  async function clickAutocompleteOption(
    component: ReactTestInstance,
    option: string
  ) {
    const autocomplete = component.find(matchesTypePredicate(Autocomplete));
    await act(async () => {
      autocomplete.props.onChange(null, option, 'select-option');
    });
  };

  function renderedValue(component: ReactTestInstance) {
    return component
      .find(matchesTypePredicate(Autocomplete)).props.value
  }

  beforeEach(() => {
    onChange = jest.fn();
    showError = jest.fn();
    input = {
      value: {picked: null},
      label: 'Abcdef',
      name: 'doesnotmatter',
      data: { 'options': ['1', '2', 'a', 'b'] },
      defaultValue: undefined,
      onChange,
      showError,
      ui_component: 'doesnotmatter',
      valueKeyMap: {picked: 'doesnotmatter'}
    };
  });

  it('mounts with label shown and using a defaultValue', () => {
    const component = mount(
      <DropdownInput
        data={input.data} value={input.value} onChange={onChange}
        label={input.label}
        defaultValue={'1'}
      />
    );

    expect(component.find(Autocomplete).length).toEqual(1);
    expect(initiallabel(component)).toEqual(input.label);
    expect(onChange).toHaveBeenCalledWith({picked: ['1']});
    expect(initialValue(component)).toEqual('1')
  });

  it('can mount without a label shown, and blank value', () => {
    const component = mount(
      <DropdownInput
        data={input.data} value={input.value} onChange={onChange}
        label={undefined}
        defaultValue={undefined}
      />
    );

    expect(component.find(Autocomplete).length).toEqual(1);
    expect(initiallabel(component)).toEqual('');
    expect(onChange).not.toHaveBeenCalled();
    expect(initialValue(component)).toEqual(null)
    expect(getInitialDropdownOptions(component)).toEqual(input.data['options'])
  });

  it('can mount without a label shown, and blank value', () => {
    const component = mount(
      <DropdownInput
        data={input.data} value={input.value} onChange={onChange}
        label={undefined}
        defaultValue={undefined}
      />
    );

    expect(component.find(Autocomplete).length).toEqual(1);
    expect(initiallabel(component)).toEqual('');
    expect(onChange).not.toHaveBeenCalled();
    expect(initialValue(component)).toEqual(null)
    expect(getInitialDropdownOptions(component)).toEqual(input.data['options'])
  });

  describe('', () => {
    beforeEach(() => {
      integrated = integrateElement(() => {
        const [value, setValue] = useState(() => input.value);
        return <DropdownInput
          onChange={(v) => {
            setValue(v);
            onChange(v);
          }}
          value={value}
          data={input.data}
          label={input.label}
          defaultValue={input.defaultValue}
        />;
      })
    });

    it('allows value selection', async () => {
      const {node} = integrated;

      expect(node.root.findAll(matchesTypePredicate(Autocomplete)).length).toEqual(1);
      expect(onChange).not.toHaveBeenCalled();

      await clickAutocompleteOption(node.root, '1')
      expect(onChange).toHaveBeenCalledWith({picked: ['1']});
      expect(renderedValue(node.root)).toEqual('1')
    });
  })
});
