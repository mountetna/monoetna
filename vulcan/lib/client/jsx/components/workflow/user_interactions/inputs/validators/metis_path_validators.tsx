import {DataEnvelope, InputValidator, ValidationInputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';
import {mapSome, maybeOfNullable, withDefault} from '../../../../../selectors/maybe';

declare const CONFIG: {[key: string]: any};

const _MetisPathValid = (
  input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>,
  type: 'file' | 'folder' | undefined,
  allow_empty_path: boolean = false
): string[] => {
  // input.value should be nested hash, like
  // {bucket: '<bucket-name>', path: 'path/to/something', type: 'file' | 'folder' | null}

  const empty_value = {bucket: '', path: '', type: null}
  const checkme = allow_empty_path ? ['bucket'] : ['bucket', 'path']

  if (!input.value || !input.value.bucket) return ['Not ready']

  const emptyKeys = withDefault(mapSome(input.value, value =>
    [...checkme].filter(o => !inputValueNonEmpty(maybeOfNullable(value[o]), true))), checkme);
  emptyKeys.sort();

  if (emptyKeys.length > 0) {
    let verb = (emptyKeys.length > 1) ? 'Values are' : 'Value is';

    return [
      `${verb} empty: ${emptyKeys.join(', ')}`
    ];
  }

  const value = withDefault(input.value, empty_value)
  // if ( (value.path as string).endsWith("/") ) {
  //   return [
  //     `Path ends in \'/\'`
  //   ];
  // }

  if (type) {
    const value = withDefault({...input.value}, empty_value)
    return type == value.type ? [] : [`Selected metis path is a ${value.type}, but a ${type} is needed.`]
  }

  return [];
};

// Is a file
export const MetisFileValidator: InputValidator<DataEnvelope<any>, DataEnvelope<any>> = (input) => {
  return _MetisPathValid(input, 'file');
};

// Is a folder
export const MetisFolderValidator: InputValidator<DataEnvelope<any>, DataEnvelope<any>> = (input) => {
  return _MetisPathValid(input, 'folder', true);
};

// Simply non-empty bucket choice
export const MetisPathValidator: InputValidator<DataEnvelope<any>, DataEnvelope<any>> = (input) => {
  return _MetisPathValid(input, undefined, true);
};
