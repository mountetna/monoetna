import React from 'react';
import {
  AllInnerValuesNotEmptyValidator,
  AllInnerValuesNotEmptyValidatorStrong
} from '../all_inner_values_not_empty_validator';
import {ValidationInputSpecification} from '../../input_types';
import {some} from "../../../../../../selectors/maybe";

describe('AllInnerValuesNotEmptyValidator', () => {
  let input: ValidationInputSpecification;

  beforeEach(() => {
    input = {
      value: null,
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
    input.value = some({});
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with null values', () => {
    input.value = some({experiment: null, tissue: null});
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with missing keys', () => {
    input.value = some({experiment: ['1']});
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with empty array values', () => {
    input.value = some({experiment: [], tissue: []});
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with array values with empty strings', () => {
    input.value = some({experiment: [''], tissue: ['']});
    expect(AllInnerValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports no errors for hashes with empty strings normally', () => {
    input.value = some({experiment: '', tissue: '1'});
    expect(AllInnerValuesNotEmptyValidator(input).length === 0).toEqual(true);
  });

    it('reports errors for hashes with empty strings when Strong', () => {
    input.value = some({experiment: '', tissue: '1'});
    expect(AllInnerValuesNotEmptyValidatorStrong(input).length > 0).toEqual(true);
  });

  it('reports no errors for well-formed hash', () => {
    input.value = some({experiment: ['1'], tissue: ['a']});
    expect(AllInnerValuesNotEmptyValidator(input).length === 0).toEqual(true);
  });
});
