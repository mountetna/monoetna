import {InputValidator, InputSpecification} from '../input_types';

const SimpleNullValidator: InputValidator = (
  input: InputSpecification
): string[] => {
  if (null == input.value) {
    return ['Value cannot be null.'];
  } else {
    return [];
  }
};

export default SimpleNullValidator;
