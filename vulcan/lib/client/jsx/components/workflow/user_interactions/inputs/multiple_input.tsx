import React, {useMemo, useCallback} from 'react';
import {DataEnvelope, InputBackendComponent, WithInputParams} from './input_types';
import InputHelp from './input_help';
import MultiselectStringInput from "./multiselect_string";
import {Maybe, maybeOfNullable, some, withDefault} from "../../../../selectors/maybe";
import {flattenNesting} from "./monoids";
import {useMemoized} from "../../../../selectors/workflow_selectors";

export function flattenOptions (
  data: {[key: string]: {[key: string]: string[]}} | null | undefined
): DataEnvelope<string[]> {
  if (data) {
    return Object.values(data).reduce((a, b) => ({...a, ...b}), {});
  }
  return {};
}

// Constructs a backend component by composing a dictionary of homogenous inner inputs, accepting each of their
// individual inner data payloads and setting a value mapping each inner input's value.
export default function MultipleInput<Value, Values>(InnerInput: InputBackendComponent<{}, Value, Values>) {
  return MultiInput;

  function MultiInput({
    value,
    data,
    onChange
  }: WithInputParams<{}, DataEnvelope<Value>, DataEnvelope<Values>>) {
    const options = useMemoized(flattenNesting, data);
    const values = useMemo(() => withDefault(value, {}), [value]);

    const onChangeValue = useCallback((label: string, v: Maybe<Value>) => {
      const newValues = {...values};
      if (v) {
        newValues[label] = v[0];
      } else {
        delete newValues[label]
      }

      onChange(some(newValues));
    }, [onChange, values])

    return (<div>
        {Object.keys(options).map(label => (<div className='view_item'>
            <div className='item_name'>{label}</div>
            <div className='item_view'>
              <InputHelp doc={label}>
                <InnerInput
                  key={label}
                  data={{options: options[label]}}
                  value={label in values ? some(values[label]) : null}
                  onChange={(value) => onChangeValue(label, value)}
                />
              </InputHelp>
            </div>
          </div>))}
      </div>);
  }
}