import {InputValidator, InputSpecification} from '../input_types';
import NotEmptyValidator from './not_empty_validator';
import {getAllOptions} from '../multiselect_string';

const ValidChoicesValidator: InputValidator = (
  input: InputSpecification
): string[] => {
  let errors = NotEmptyValidator(input);

  const options = getAllOptions(input.data);

  if (
    input &&
    input.value &&
    !Object.values(input.value).every((v) => options.includes(v as string))
  ) {
    errors.push('Invalid choice selected.');
  }

  return errors;
};

export default ValidChoicesValidator;
