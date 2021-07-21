// Input component that takes a simple object and
// shows values based on a selected key
import React, {useCallback, useMemo, useState,} from 'react';

import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import CheckboxesInput from './checkboxes';
import {DataEnvelope, WithInputParams} from './input_types';
import {Maybe, some} from "../../../../selectors/maybe";
import {joinNesting} from "./monoids";
import {useMemoized} from "../../../../selectors/workflow_selectors";
import {useSetsDefault} from "./useSetsDefault";

function DropdownCheckboxCombo({
  dropdownValue,
  handleSelect,
  checkboxOptions,
  dropdownOptions,
  value,
  onChange,
}: {
  dropdownValue: string | null;
  checkboxOptions: string[];
  dropdownOptions: string[];
  handleSelect: (value: string) => void;
  value: string[];
  onChange: (values: Maybe<string[]>) => void;
}) {
  return (
    <React.Fragment>
      <div>
        <DropdownAutocomplete
          value={dropdownValue}
          onSelect={handleSelect}
          list={dropdownOptions}
        />
      </div>
      <div className='checkbox-input-wrapper'>
        <CheckboxesInput
          data={{checkboxOptions}}
          onChange={onChange}
          value={some(value)}
        />
      </div>
    </React.Fragment>
  );
}

export default function SingleDropdownMulticheckbox({ data, onChange, ...props }: WithInputParams<{}, DataEnvelope<string[]>, DataEnvelope<string[]>>) {
  // input.data for this component should be
  // {'a': {experiment: ['1', '2'], tissue: ['a', 'b']}}
  const allOptions = useMemoized(joinNesting, data);
  const [dropdownValue, setDropdownValue] = useState(() => Object.keys(allOptions)[0] as string | null);
  const value = useSetsDefault(allOptions, props.value, onChange);
  const dropdownOptions = useMemo(() => Object.keys(allOptions), [allOptions]);
  const checkboxOptions = dropdownValue == null ? [] : allOptions[dropdownValue];
  const selectedValue = dropdownValue == null ? [] : value[dropdownValue] || [];
  const onChangeInner = useCallback((values: Maybe<string[]>) => {
    if (!dropdownValue || values == null) return;
    onChange(some({
      ...value,
      [dropdownValue]: [...values[0]],
    }))
  }, [dropdownValue, onChange, value]);

  return (
    <div>
      <DropdownCheckboxCombo
        checkboxOptions={checkboxOptions}
        dropdownOptions={dropdownOptions}
        dropdownValue={dropdownValue}
        handleSelect={setDropdownValue}
        value={selectedValue}
        onChange={onChangeInner}
      />
    </div>
  );
};
