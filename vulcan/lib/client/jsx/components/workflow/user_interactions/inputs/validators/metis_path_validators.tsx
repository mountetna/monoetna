import {DataEnvelope, InputValidator, ValidationInputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';
import {mapSome, maybeOfNullable, withDefault} from '../../../../../selectors/maybe';

declare const CONFIG: {[key: string]: any};

const _MetisPathValid = (
  input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>,
  type: 'file' | 'folder' | undefined,
  allow_toplevel_path: boolean = false
): string[] => {
  // input.value should be nested hash, like
  // {bucket: '<bucket-name>', path: 'path/to/something', type: 'file' | 'folder' | null}

  const empty_value = {bucket: '', path: '', type: null}
  const checkme = allow_toplevel_path ? ['bucket'] : ['bucket', 'path']

  const value = withDefault(input.value, undefined)

  if (!value) return ['Value not ready']

  // Missing bucket or path?
  const emptyKeys = withDefault(mapSome(input.value, value =>
    [...checkme].filter(o => !inputValueNonEmpty(maybeOfNullable(value[o]), true))), checkme);
  emptyKeys.sort();

  if (emptyKeys.length > 0) {
    let verb = (emptyKeys.length > 1) ? 'Values are' : 'Value is';

    return [
      `${verb} empty: ${emptyKeys.join(', ')}`
    ];
  }

  if (type) {
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
