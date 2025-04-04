import {DataEnvelope, InputValidator, ValidationInputSpecification} from '../../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';
import {mapSome, maybeOfNullable, withDefault} from '../../../../../selectors/maybe';

declare const CONFIG: {[key: string]: any};

const _MetisPathValid = (
  input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>,
  type: 'file' | 'folder' | undefined,
  allow_toplevel_path: boolean = false,
  path_regex_string?: string,
  path_regex_descriptor?: string
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

  if (type && type != value.type) {
    return [`Selected metis path is a ${value.type}, but a ${type} is needed.`]
  }

  if (path_regex_string) {
    const path_regex = new RegExp(path_regex_string)
    if (!path_regex.test(value.path)) {
      const description = path_regex_descriptor!=null ?
        `should target a ${path_regex_descriptor} (by matching the regex ${path_regex})` :
        `must match the regex ${path_regex}`
      return [
        `Selected metis path ${description}`
      ];
    }
  }

  return [];
};

// Is a file
export function MetisFileValidator(path_regex_string?: string, path_regex_descriptor?: string): InputValidator<DataEnvelope<any>, DataEnvelope<any>> {
  return (input) => {
    return _MetisPathValid(input, 'file', false, path_regex_string, path_regex_descriptor);
  };
}

// Is a folder
export function MetisFolderValidator(path_regex_string?: string, path_regex_descriptor?: string): InputValidator<DataEnvelope<any>, DataEnvelope<any>> {
  return (input) => {
    return _MetisPathValid(input, 'folder', true, path_regex_string, path_regex_descriptor);
  };
}

// Simply non-empty bucket choice
export function MetisPathValidator(path_regex_string?: string, path_regex_descriptor?: string): InputValidator<DataEnvelope<any>, DataEnvelope<any>> {
  return (input) => {
    return _MetisPathValid(input, undefined, true, path_regex_string, path_regex_descriptor);
  };
}
