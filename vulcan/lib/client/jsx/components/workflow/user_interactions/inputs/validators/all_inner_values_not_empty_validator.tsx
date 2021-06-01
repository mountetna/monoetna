import * as _ from 'lodash';

import {InputValidator, InputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';
import {flattenOptions} from '../multiple_multiselect_string_all';

const AllInnerValuesNotEmptyValidator: InputValidator = (
  input: InputSpecification
): string[] => {
  // input.value should be nested hash, like
  // {experiment: [list,of,options], tissue: [other,choices]}
  const msg = 'All inner values must be selected.';
  if (null == input.value || Object.values(input.value).length === 0)
    return [msg];

  if (
    !_.isEqual(
      Object.keys(input.value),
      Object.keys(flattenOptions(input.data))
    )
  ) {
    return [msg];
  } else if (
    !Object.values(input.value).every((selections) =>
      inputValueNonEmpty(selections)
    )
  ) {
    return [msg];
  } else {
    return [];
  }
};

export default AllInnerValuesNotEmptyValidator;
