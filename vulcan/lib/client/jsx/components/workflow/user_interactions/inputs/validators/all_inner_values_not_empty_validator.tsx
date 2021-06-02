import * as _ from 'lodash';

import {InputValidator, InputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';
import {flattenOptions} from '../multiple_multiselect_string_all';
import NotEmptyValidator from './not_empty_validator';

const AllInnerValuesNotEmptyValidator: InputValidator = (
  input: InputSpecification
): string[] => {
  // input.value should be nested hash, like
  // {experiment: [list,of,options], tissue: [other,choices]}
  if (null == input.value || Object.values(input.value).length === 0)
    return ['All inner values must be selected.'];

  if (
    !_.isEqual(
      Object.keys(input.value).sort(),
      Object.keys(flattenOptions(input.data)).sort()
    )
  ) {
    console.error(
      Object.keys(input.value).sort(),
      Object.keys(flattenOptions(input.data)).sort()
    );
    return ['Missing at least one inner value key.'];
  } else if (
    !Object.values(input.value).every((selections) =>
      inputValueNonEmpty(selections)
    )
  ) {
    console.error(Object.values(input.value));
    return ['At least one inner value key is empty.'];
  }

  return [];
};

export default AllInnerValuesNotEmptyValidator;
