import React, {useCallback} from 'react';
import {DataEnvelope, InputBackendComponent, WithInputParams} from './input_types';
import InputHelp from './input_help';
import {Maybe, some} from "../../../../selectors/maybe";
import {joinNesting} from "./monoids";
import {useMemoized} from "../../../../selectors/workflow_selectors";
import {useSetsDefault} from "./useSetsDefault";

// Constructs a backend component by composing a dictionary of homogenous inner inputs, accepting each of their
// individual inner data payloads and setting a value mapping each inner input's value.
export default function MultipleInput<Value, Values>(InnerInput: InputBackendComponent<{}, Value, Values>) {
  return MultiInput;

  function MultiInput({
    data,
    onChange,
    ...props
  }: WithInputParams<{}, DataEnvelope<Value>, DataEnvelope<Values>>) {
    const options = useMemoized(joinNesting, data);
    const values = useSetsDefault({}, props.value, onChange);

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