import {InputValidator, InputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';

const AllInnerValuesNotEmptyValidator: InputValidator = (
  input: InputSpecification
): string[] => {
  const msg = 'All inner values must be selected.';
  if (null == input.value || Object.values(input.value).length === 0)
    return [msg];

  if (
    // input.value should be nested hash, like
    // {experiment: [list,of,options], tissue: [other,choices]}
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
