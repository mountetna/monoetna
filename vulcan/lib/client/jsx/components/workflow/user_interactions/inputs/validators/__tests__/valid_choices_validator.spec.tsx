import React from 'react';
import ValidChoicesValidator from '../valid_choices_validator';
import {InputSpecification} from '../../input_types';

describe('ValidChoicesValidator', () => {
  let input: InputSpecification;

  beforeEach(() => {
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
      name: 'test-input',
      data: ['1', '2', 'a', 'b']
    };
  });

  it('reports errors for null values', () => {
    expect(ValidChoicesValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for invalid choice selection', () => {
    input.value = ['x'];
    expect(ValidChoicesValidator(input).length > 0).toEqual(true);
  });

  it('reports no errors for valid choice selection', () => {
    input.value = ['a'];
    expect(ValidChoicesValidator(input).length === 0).toEqual(true);
  });
});
