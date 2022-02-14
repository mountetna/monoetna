import {
  DataEnvelope,
  InputValidator,
  ValidationInputSpecification
} from '../input_types';
import {isSome, withDefault} from '../../../../../selectors/maybe';
import {isNullish} from '../../../../../selectors/workflow_selectors';

const _AllInnerKeysNotNullValidator = (
  input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>,
  innerKeyName = 'data'
): string[] => {
  // input.value should be nested hash, like
  // {data: {[key: string]: any}, sourceData: {[key: string]: any}}
  // This validator checks that no key == `null` within an inner Hash.
  if (!isSome(input.value)) return ['Input value is null.'];

  const innerHash = withDefault(input.value, {});

  if (!innerHash.hasOwnProperty(innerKeyName))
    return [`Inner hash missing key ${innerKeyName}.`];

  const nullKeys = Object.keys(innerHash[innerKeyName]).filter(isNullish);

  if (nullKeys.length > 0) {
    let noun = nullKeys.length > 1 ? 'keys' : 'key';

    return [`Found ${nullKeys.length} null ${noun}.`];
  }

  return [];
};

export const AllInnerKeysNotNullValidator: InputValidator<
  DataEnvelope<any>,
  DataEnvelope<any>
> = (input) => {
  return _AllInnerKeysNotNullValidator(input);
};

export default AllInnerKeysNotNullValidator;
