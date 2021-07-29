import React from 'react';
import NotEmptyValidator from '../not_empty_validator';
import {ValidationInputSpecification} from '../../input_types';
import {some} from "../../../../../../selectors/maybe";

describe('NotEmptyValidator', () => {
  let input: ValidationInputSpecification<any>;

  beforeEach(() => {
    input = {
      value: null,
      data: {},
    };
  });

  it('reports errors for unset values', () => {
    expect(NotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for null values', () => {
    input.value = some(null);
    expect(NotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for empty arrays', () => {
    input.value = some([]);
    expect(NotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for arrays only containing an empty string', () => {
    input.value = some(['']);
    expect(NotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports no errors for populated arrays', () => {
    input.value = some(['choice 1']);
    expect(NotEmptyValidator(input).length === 0).toEqual(true);
  });

  it('reports no errors for non-null values', () => {
    input.value = some('something!');
    expect(NotEmptyValidator(input).length === 0).toEqual(true);
  });
});
