import React from 'react';
import NotEmptyValidator from '../not_empty_validator';
import {InputSpecification} from '../../input_types';

describe('NotEmptyValidator', () => {
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
    expect(NotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for empty arrays', () => {
    input.value = [];
    expect(NotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for arrays with empty string only', () => {
    input.value = [''];
    expect(NotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports no errors for populated arrays', () => {
    input.value = ['choice 1'];
    expect(NotEmptyValidator(input).length === 0).toEqual(true);
  });

  it('reports no errors for non-null values', () => {
    input.value = 'something!';
    expect(NotEmptyValidator(input).length === 0).toEqual(true);
  });
});
