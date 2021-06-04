import React from 'react';
import AllInnerValuesNotEmptyValidator from '../all_inner_values_not_empty_validator';
import {InputSpecification} from '../../input_types';

describe('AllInnerValuesNotEmptyValidator', () => {
  let input: InputSpecification;

  beforeEach(() => {
    input = {
      type: 'doesnotmatter',
      value: null,
      label: 'Abcdef',
      name: 'test-input',
      data: {
        a: {
          experiment: ['1', '2'],
          tissue: ['a', 'b']
        }
      }
    };
  });

  it('reports errors for null values', () => {
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for empty hashes', () => {
    input.value = {};
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with null values', () => {
    input.value = {experiment: null, tissue: null};
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with missing keys', () => {
    input.value = {experiment: ['1']};
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with empty array values', () => {
    input.value = {experiment: [], tissue: []};
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with array values with empty strings', () => {
    input.value = {experiment: [''], tissue: ['']};
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports no errors for well-formed hash', () => {
    input.value = {experiment: ['1'], tissue: ['a']};
    expect(AllInnerValuesNotEmptyValidator(input).length === 0).toEqual(true);
  });
});
