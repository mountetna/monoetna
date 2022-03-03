import {
  DataEnvelope,
  InputValidator,
  ValidationInputSpecification
} from '../input_types';
import {isSome, withDefault} from '../../../../../selectors/maybe';
import {isNullish} from '../../../../../selectors/workflow_selectors';

const _AllInnerKeysNotNullValidator = (
  input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>
): string[] => {
  // input.value should be nested hash, like
  // {data: {[key: string]: any}, sourceData: {[key: string]: any}}
  // This validator checks that no key == `null` within an inner Hash.
  if (!isSome(input.value)) return ['Input value is null.'];

  const innerHash = withDefault(input.value, {});
  const innerKeys = Object.keys(innerHash);

  if (innerKeys.length === 0) return ['Inner hash is empty!'];

  if (innerKeys.some((key) => isNullish(innerHash[key])))
    return ['Contains undefined inner hashes.'];

  let errors: string[] = [];

  innerKeys.forEach((key) => {
    const nullKeys = Object.keys(innerHash[key]).filter(isNullish);

    if (nullKeys.length > 0) {
      let noun = nullKeys.length > 1 ? 'keys' : 'key';

      errors.push(`Found ${nullKeys.length} null ${noun}.`);
    }
  });
  return errors;
};

export const AllInnerKeysNotNullValidator: InputValidator<
  DataEnvelope<any>,
  DataEnvelope<any>
> = (input) => {
  return _AllInnerKeysNotNullValidator(input);
};

export default AllInnerKeysNotNullValidator;
