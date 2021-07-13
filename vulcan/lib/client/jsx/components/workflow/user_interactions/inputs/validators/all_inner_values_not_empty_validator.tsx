import {InputValidator, ValidationInputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';
import {flattenOptions} from '../multiple_multiselect_string_all';
import {mapSome, maybeOfNullable, withDefault} from "../../../../../selectors/maybe";

const _AllInnerValuesNotEmptyValidator = (
  input: ValidationInputSpecification, strong = false
): string[] => {
  // input.value should be nested hash, like
  // {experiment: [list,of,options], tissue: [other,choices]}
  const allOptions = flattenOptions(input.data);

  const emptyKeys = withDefault(mapSome(input.value, value =>
    Object.keys(allOptions).filter(o => inputValueNonEmpty(maybeOfNullable(o), strong))), Object.keys(allOptions));
  emptyKeys.sort();

  if (emptyKeys.length > 0) {
    let verb = (emptyKeys.length > 1) ? 'Values are' : 'Value is';

    return [
      `${verb} empty: ${emptyKeys.join(',')}`
    ];
  }

  return [];
};

export const AllInnerValuesNotEmptyValidator: InputValidator = (input) => {
  return _AllInnerValuesNotEmptyValidator(input, false);
};

export const AllInnerValuesNotEmptyValidatorStrong: InputValidator = (input) => {
  return _AllInnerValuesNotEmptyValidator(input, true);
};

export default AllInnerValuesNotEmptyValidator;
