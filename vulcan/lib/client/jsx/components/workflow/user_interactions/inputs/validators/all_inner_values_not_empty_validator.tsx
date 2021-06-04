import * as _ from 'lodash';

import {InputValidator, InputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';
import {flattenOptions} from '../multiple_multiselect_string_all';

const AllInnerValuesNotEmptyValidator: InputValidator = (
  input: InputSpecification
): string[] => {
  // input.value should be nested hash, like
  // {experiment: [list,of,options], tissue: [other,choices]}
  const allOptions = flattenOptions(input.data);

  function findEmptyKeys() {
    if (null == input.value || Object.values(input.value).length === 0) {
      return Object.keys(allOptions).sort();
    }

    return Object.keys(allOptions)
      .filter(
        (o: string) => !(o in input.value) || !inputValueNonEmpty(input.value[o])
      ).sort()
  }

  let emptyKeys = findEmptyKeys();

  if (emptyKeys.length > 0) {
    let verb = (emptyKeys.length > 1) ? 'Values are' : 'Value is';

    return [
      `${verb} empty: ${emptyKeys.join(',')}`
    ];
  }

  return [];
};

export default AllInnerValuesNotEmptyValidator;
