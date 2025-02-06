import {DataEnvelope, InputValidator, ValidationInputSpecification} from '../../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';
import {mapSome, maybeOfNullable, withDefault} from '../../../../../selectors/maybe';
import {joinNesting} from '../../monoids';

const _AllInnerValuesNotEmptyValidator = (
  input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>, disallow_empty_string = false, disallow_empty_array = true
): string[] => {
  // input.value should be nested hash, like
  // {experiment: [list,of,options], tissue: [other,choices]}
  const allOptions = joinNesting(input.data);

  const emptyKeys = withDefault(mapSome(input.value, value =>
    Object.keys(allOptions).filter(o => !inputValueNonEmpty(maybeOfNullable(value[o]), disallow_empty_string, disallow_empty_array))), Object.keys(allOptions));
  emptyKeys.sort();

  if (emptyKeys.length > 0) {
    let verb = (emptyKeys.length > 1) ? 'Values are' : 'Value is';

    return [
      `${verb} empty: ${emptyKeys.join(', ')}`
    ];
  }

  return [];
};

export const AllInnerValuesNotEmptyValidator: InputValidator<DataEnvelope<any>, DataEnvelope<any>> = (input) => {
  return _AllInnerValuesNotEmptyValidator(input, false);
};

export const AllInnerValuesNotEmptyValidatorStrong: InputValidator<DataEnvelope<any>, DataEnvelope<any>> = (input) => {
  return _AllInnerValuesNotEmptyValidator(input, true);
};

export const AllInnerValuesNotEmptyAllowingEmptyArrayValidator: InputValidator<DataEnvelope<any>, DataEnvelope<any>> = (input) => {
  return _AllInnerValuesNotEmptyValidator(input, false, false);
};

export const AllInnerValuesNotEmptyAllowingEmptyArrayValidatorStrong: InputValidator<DataEnvelope<any>, DataEnvelope<any>> = (input) => {
  return _AllInnerValuesNotEmptyValidator(input, true, false);
};

export default AllInnerValuesNotEmptyValidator;
