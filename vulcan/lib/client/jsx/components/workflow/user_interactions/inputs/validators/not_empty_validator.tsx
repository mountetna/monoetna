import {InputValidator, InputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';

const NotEmptyValidator: InputValidator = (
  input: InputSpecification
): string[] => {
  if (!inputValueNonEmpty(input.value)) {
    return ['Value cannot be empty.'];
  } else {
    return [];
  }
};

export default NotEmptyValidator;
