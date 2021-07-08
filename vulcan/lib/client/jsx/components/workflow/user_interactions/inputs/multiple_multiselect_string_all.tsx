import React, {useMemo, useCallback} from 'react';
import {WithInputParams} from './input_types';
import InputHelp from './input_help';
import MultiselectStringInput from "./multiselect_string";
import {Maybe, maybeOfNullable, some, withDefault} from "../../../../selectors/maybe";

export const flattenOptions = (
  data: {[key: string]: {[key: string]: string[]}} | null | undefined
) => {
  if (data) {
    return Object.values(data).reduce((a, b) => ({...a, ...b}), {});
  }
  return {};
};

export default function MultipleMultiselectStringAllInput({
  value,
  data,
  onChange
}: WithInputParams<{}>) {
  const options: {[k: string]: string[]} = useMemo(() => flattenOptions(data), [
    data
  ]);
  const values: {[k: string]: string} = useMemo(() => withDefault(value, {}), [value]);
  const onChangeValue = useCallback((label: string, v: Maybe<string>) => {
    // this does not work, come back to it.
    const newValues: {[k: string]: string} = {...values};
    if (v) {
      newValues[label] = v[0];
    } else {
      delete newValues[label]
    }

    onChange(some(newValues));
  }, [onChange, values])

  return (
    <div>
      {Object.keys(options).map(label => (
        <div className='view_item'>
          <div className='item_name'>{label}</div>
          <div className='item_view'>
            <InputHelp doc={label}>
              <MultiselectStringInput
                key={label}
                data={{options}}
                value={maybeOfNullable(values[label])}
                onChange={(value) => onChangeValue(label, value)}
              />
            </InputHelp>
          </div>
        </div>
      ))}
    </div>
  );
};
