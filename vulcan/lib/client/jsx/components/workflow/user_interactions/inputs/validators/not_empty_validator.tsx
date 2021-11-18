import {InputValidator, BoundInputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';

export const NotEmptyValidator: InputValidator = (input) => {
  if (!inputValueNonEmpty(input.value)) {
    return ['Value cannot be empty.'];
  } else {
    return [];
  }
};

export const StronglyNotEmptyValidator: InputValidator = (input) => {
  if (!inputValueNonEmpty(input.value, true)) {
    return ['Value cannot be empty.'];
  } else {
    return [];
  }
};
