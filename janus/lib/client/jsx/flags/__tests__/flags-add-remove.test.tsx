import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import AddRemoveFlag from '../flags-add-remove';

import {FlagsState} from '../flags-context';

describe('AddRemoveFlag', () => {
  let state: FlagsState;

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
  });

  it('renders correctly', () => {
    const component = mount(<FlagsView />);
  });

  it('shows Add/Remove card when users are selected', () => {});

  it('filters users based on text search', () => {});

  it('filters users based on project search', () => {});

  it('filters users based on flag search', () => {});

  it('filters users based on all 3 search params', () => {});
});
