import AllOutputValuesNotEmptyValidator from '../all_output_values_not_empty_validator';
import {DataEnvelope, ValidationInputSpecification} from '../../input_types';
import {some} from "../../../../../../selectors/maybe";

describe('AllOutputValuesNotEmptyValidator', () => {
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
    expect(AllOutputValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for empty hashes', () => {
    input.value = some({});
    expect(AllOutputValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with null values', () => {
    input.value = some({experiment: null, tissue: null});
    expect(AllOutputValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports no errors for hashes with missing keys', () => {
    input.value = some({experiment: ['1']});
    expect(AllOutputValuesNotEmptyValidator(input)).toEqual([]);
  });

  it('reports errors for hashes with empty array values', () => {
    input.value = some({experiment: [], tissue: []});
    expect(AllOutputValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports errors for hashes with array values with empty strings', () => {
    input.value = some({experiment: [''], tissue: ['']});
    expect(AllOutputValuesNotEmptyValidator(input).length > 0).toEqual(true);
  });

  it('reports no errors for hashes with empty strings normally', () => {
    input.value = some({experiment: '', tissue: '1'});
    expect(AllOutputValuesNotEmptyValidator(input)).toEqual([]);
  });

  it('reports no errors for well-formed hash', () => {
    input.value = some({experiment: ['1'], tissue: ['a']});
    expect(AllOutputValuesNotEmptyValidator(input).length === 0).toEqual(true);
  });
});
