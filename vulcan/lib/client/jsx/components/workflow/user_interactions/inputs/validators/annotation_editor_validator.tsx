import {
  DataEnvelope,
  InputValidator,
  ValidationInputSpecification
} from '../input_types';
import AllInnerKeysNotNullValidator from './all_inner_keys_not_null_validator';
import { joinNesting } from '../monoids';

const _AnnotationEditorValidator = (
  input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>
): string[] => {

  // First, validating that all elements are filled
  let errors = AllInnerKeysNotNullValidator(input);
  if (errors.length>0) return errors;

  // Check column names
  //   From above, already know input.value is filled, thus withDefault not needed
  //   The 2 outputs, formulaic_data and calculated_data, can be in either order, but both will have the same keys & key order, so we don't care which we get
  const outputDF = Object.values((input.value as object[])[0])[0] as object | object[];
  // For some reason, the output structure changes when there is one column versus multiple
  const userColNames = Array.isArray(outputDF) ? Object.keys(outputDF[0]) : Object.keys(outputDF);
  const inputDF = joinNesting(input.data);

  // Ensuring the clustering column name (first column) matches its original value
  if (Object.keys(inputDF)[0] != userColNames[0]) {
    errors.push('Clustering metadata name has been modified. This would break annotation import steps.')
  };

  // Ensure at least one column starts with 'annot'
  if (userColNames.filter((val) => val.startsWith('annot')).length < 1) {
    errors.push('No annotation columns (named starting with \'annot\') detected. At least 1 is required.')
  };

  return errors;
};

export const AnnotationEditorValidator: InputValidator<
  DataEnvelope<any>,
  DataEnvelope<any>
> = (input) => {
  return _AnnotationEditorValidator(input);
};

export default AnnotationEditorValidator;
