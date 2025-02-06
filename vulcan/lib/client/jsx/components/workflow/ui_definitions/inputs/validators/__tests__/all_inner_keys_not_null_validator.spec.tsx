import React from 'react';
import {AllInnerKeysNotNullValidator} from '../all_inner_keys_not_null_validator';
import {DataEnvelope, ValidationInputSpecification} from '../../input_types';
import {some} from '../../../../../../selectors/maybe';

describe('AllInnerKeysNotNullValidator', () => {
  let input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>;

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
    expect(AllInnerKeysNotNullValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for empty hashes', () => {
    input.value = some({});
    expect(AllInnerKeysNotNullValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with null values', () => {
    input.value = some({
      data: {
        null: {}
      }
    });
    expect(AllInnerKeysNotNullValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with undefined values', () => {
    input.value = some({
      data: undefined
    });
    expect(AllInnerKeysNotNullValidator(input).length > 0).toEqual(true);
  });

  it('reports no errors for well-formed hash', () => {
    input.value = some({data: {column1: {'1': 1, '0': 0}}});
    expect(AllInnerKeysNotNullValidator(input).length === 0).toEqual(true);
  });
});
