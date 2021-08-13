import { inputValueNonEmpty } from '../../../../../selectors/workflow_selectors';
import {DataEnvelope, ValidationInputSpecification} from '../input_types';
import AllInnerValuesNotEmptyValidator, { AllInnerValuesNotEmptyValidatorStrong } from './all_inner_values_not_empty_validator';
import NotEmptyValidator from './not_empty_validator';

const AllOutputValuesNotEmptyValidator = (input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>): string[] => {
  // When this is used, we don't care if keys in input.data are missing from input.value
  if (!inputValueNonEmpty(input.value)) {
    return ['Value cannot be empty.'];
  } else {
    return AllInnerValuesNotEmptyValidator({...input, data: input.value as DataEnvelope<any>});
  }
};

export default AllOutputValuesNotEmptyValidator;