import { inputValueNonEmpty } from '../../../../../selectors/workflow_selectors';
import {DataEnvelope, ValidationInputSpecification} from '../input_types';
import AllInnerValuesNotEmptyValidator, { AllInnerValuesNotEmptyAllowingEmptyArrayValidator } from './all_inner_values_not_empty_validator';

const AllOutputValuesNotEmptyValidator = (input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>): string[] => {
  // When this is used, we don't care if keys in input.data are missing from input.value
  if (!inputValueNonEmpty(input.value)) {
    return ['Value cannot be empty.'];
  } else {
    return AllInnerValuesNotEmptyValidator({...input, data: input.value as DataEnvelope<any>});
  }
};

export const AllOutputValuesNotEmptyAllowingEmptyArrayValidator = (input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>): string[] => {
  // Compatively, allows empty arrays.
  if (!inputValueNonEmpty(input.value)) {
    return ['Value cannot be empty.'];
  } else {
    return AllInnerValuesNotEmptyAllowingEmptyArrayValidator({...input, data: input.value as DataEnvelope<any>});
  }
};

export default AllOutputValuesNotEmptyValidator;