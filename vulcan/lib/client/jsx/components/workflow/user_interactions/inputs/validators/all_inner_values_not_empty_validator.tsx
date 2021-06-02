import * as _ from 'lodash';

import {InputValidator, InputSpecification} from '../input_types';
import {inputValueNonEmpty} from '../../../../../selectors/workflow_selectors';
import {flattenOptions} from '../multiple_multiselect_string_all';

const AllInnerValuesNotEmptyValidator: InputValidator = (
  input: InputSpecification
): string[] => {
  // input.value should be nested hash, like
  // {experiment: [list,of,options], tissue: [other,choices]}
  if (null == input.value || Object.values(input.value).length === 0)
    return ['All inner values must be selected.'];

  const allOptions = flattenOptions(input.data);

  function findInvalidSelections() {
    return Object.entries(input.value).filter(
      ([key, selections]: [string, any]) => {
        let validOptions = allOptions[key];
        if (!validOptions) return true;

        return !selections.every((s: string) => validOptions.includes(s));
      }
    );
  }

  function findEmptyKeys() {
    return Object.values(input.value).filter(
      (selections) => !inputValueNonEmpty(selections)
    );
  }

  function findMissingKeys() {
    if (null == input.value) return Object.keys(allOptions).sort();

    return Object.keys(allOptions)
      .filter((o) => !Object.keys(input.value).includes(o))
      .sort();
  }

  if (findMissingKeys().length > 0) {
    console.error(
      Object.keys(input.value).sort(),
      Object.keys(allOptions).sort()
    );
    let missingKeys = findMissingKeys();
    return [
      `Missing inner value key${
        missingKeys.length > 1 ? 's' : ''
      }: ${missingKeys.map((k) => k).join(',')}`
    ];
  } else if (findEmptyKeys().length > 0) {
    console.error(Object.values(input.value));
    let emptyKeys = findEmptyKeys();
    let verb = 'value is';
    if (emptyKeys.length > 1) {
      verb = 'values are';
    }

    return [
      `Inner key ${verb} empty: ${emptyKeys.map((k: any) => k).join(',')}`
    ];
  } else if (findInvalidSelections().length > 0) {
    let invalidSelections = findInvalidSelections();
    return [
      `Invalid selection${
        invalidSelections.length > 1 ? 's' : ''
      } present: ${invalidSelections
        .map(
          ([key, selections]: [string, any]) =>
            `${key} -- ${selections.join(',')}`
        )
        .join(',')}`
    ];
  }

  return [];
};

export default AllInnerValuesNotEmptyValidator;
