// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useMemo, useState} from 'react';
import * as _ from 'lodash';

import {DataEnvelope, WithInputParams} from './input_types';
import { useSetsDefault } from './useSetsDefault';
import { maybeOfNullable, some, withDefault, Maybe } from '../../../../selectors/maybe';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import MultiselectStringInput from './multiselect_string';
import { Button, Slider } from '@material-ui/core';
import { pick } from 'lodash';
import StringInput from './string';
import BooleanInput from './boolean';
import { flattenStringOptions } from './monoids';

/*
This input is closely tied to a python script for running scanpy-based differential expression.

Input:
- The metadata data.frame, as json, of an single-cell object.

Output:
- A dictionary containing:
    - method = a string of "btwn-all-de-groups", "btwn-groups", "btwn-groups-multiple"
    - (optional) subset_meta = a string naming a metadata to subset by (when: subsetting needed)
    - (optional) subset_use = a string vector of values of the subset meta to keep (when: subsetting needed)
    - de_meta = a string naming a metadata to compare between
    - de_group_1 = a string vector of values of the de_meta to either compare among (when: "btwn-all-de-groups") or use as group1 (when:"btwn-groups", "btwn-groups-multiple")
    - (optional) de_group_2 = a string vector of values of the de_meta to use as group2 (when: "btwn-groups", "btwn-groups-multiple")
    - (optional) groups_meta = a string naming a metadata to use as cell groupings (when: "btwn-groups-multiple")
    - (optional) groups_use = a string vector indicating which values of the groupings meta to investigate within (when: "btwn-groups-multiple")

Labels:
- Between all labels of a particular metadata
- Between particular labels
- Between particular labels in multiple cell groups
*/

export default function DiffExpSC({
  data, onChange, ...props
}: WithInputParams<{}, DataEnvelope<any>, any>, plotType: string) {
  const value = useSetsDefault({method: null}, props.value, onChange);

  const options: string[] = useMemo(() => {
    if (data == null) return [];
    return Object.keys(data)
    // ToDo: only return keys where values are strings
  }, [data]);
  
  const subOptions: DataEnvelope<string[]> = useMemo(() => {
    if (data == null) return [];
    let opts = {...data}
    Object.entries(opts).map( ([key,value]) => {
      opts[key] = flattenStringOptions({'0': value})
    })
    return opts;
  }, [data])
  
  const addSubset = ()
  
  const removeSubset = (vals = {...props.value} as DataEnvelope<any>) => {
    if (Object.keys(vals).includes('subset_meta')) delete vals['subset_meta']
    if (Object.keys(vals).includes('subset_use')) delete vals['subset_use']
    onChange(vals);
  }

  const updateValue = (newValue: any, key: string, prevValues = {...value}) => {
    prevValues[key] = newValue;
    onChange(some(prevValues));
  };

  /*
  method question
  subset button
    subset_by
    subset_use
  DE_by
    DE_group1
    DE_group2
  group_by if there
  */
  
  return (
    <div>
      {Object.entries(shownValues).map(([key, val]) => {
        return component_use(key, val, extra_inputs[key])
      })}
      <div>
        <Button
          variant="contained"
          color="primary"
          onClick={() => {toggleAdvanced()}}
          >
          {showHide} Advanced Options
        </Button>
      </div>
    </div>
  );

};

const method_labels: DataEnvelope<string> = {
  'btwn-all-de-groups': 'Between all labels of a particular metadata',
  'btwn-groups': 'Between particular labels',
  'btwn-groups-multiple': 'Between particular labels in multiple cell groups'
}

const output_sets: DataEnvelope<DataEnvelope<string|string[]|null>> = {
  'btwn-all-de-groups': {
    'method': 'btwn-all-de-groups',
    'de_meta': null,
    'de_group_1': null
  },
  'btwn-groups': {
    'method': 'btwn-all-de-groups',
    'de_meta': null,
    'de_group_1': null,
    'de_group_2': null
  },
  'btwn-groups-multiple': {
    'method': 'btwn-all-de-groups',
    'de_meta': null,
    'de_group_1': null,
    'de_group_2': null,
    'groups_meta': null,
    'groups_use': null
  }
}

// Component Setups
const stringInput = (
  key: string = "filler", changeFxn: Function, value: string = "filler",
  label: string = 'hello') => {
    return (
      <StringInput
        key={key}
        label={label}
        data={val_wrap(value)}
        value={maybeOfNullable(value)}
        onChange={(newValue) => changeFxn(withDefault(newValue,'make'), key)}
      />
    )
  };

// const nestableDropdownInput = (
//   key: string = "filler", changeFxn: Function, value: string | null,
//   label: string, options: DataEnvelope<null>) => {
//    
//     return(
//       <NestedSelectAutocompleteInput
//         key={key}
//         label={label}
//         data={{options}} 
//         value={maybeOfNullable(value)}
//         onChange={(val) => changeFxn(withDefault(val, null), key)}
//       />
//     )
//   }

const checkboxInput = (
  key: string = "filler", changeFxn: Function, value: boolean = false,
  label: string) => {

    return(
      <BooleanInput
        key={key}
        label={label}
        value={maybeOfNullable(value)}
        data={val_wrap(value)}
        onChange={ value => changeFxn(withDefault(value,false), key)}
      />
    )
  }

const dropdownInput = (
  key: string = "filler", changeFxn: Function, value: string | null,
  label: string, options: string[]) => {
    
    return(
      <div key={key}>
        {label}
        <DropdownAutocomplete
          list={options}
          value={value}
          onSelect={(val: string) => changeFxn(val, key)}
        />
      </div>
    )
  }

const MultiselectInput = (
  key: string = "filler", changeFxn: Function, value: string[],
  label: string, options: string[]) => {
    
    return(
      <div key={key}>
        {label}
        <MultiselectStringInput
          data={{'0': options}}
          value={some(value)}
          onChange={(val: Maybe<string[]>) => changeFxn(withDefault(val, null), key)}
        />
      </div>
    )
  }

const sliderInput = (
  key: string = "filler", changeFxn: Function, value: number,
  label: string, min: number = 0.1, max: number = 20) => {

    return(
        <div key={key}>
          {label}
          <Slider
            value={value}
            onChange={(event, newValue) => changeFxn(newValue as number, key)}
            min={min}
            max={max}
            valueLabelDisplay="auto"
          />
        </div>
    )
  }

const comps: DataEnvelope<Function> = {
  'plot_title': stringInput,
  'legend_title': stringInput,
  'xlab': stringInput,
  'ylab': stringInput,
  'x_by': dropdownInput,
  'y_by': dropdownInput,
  'color_by': dropdownInput,
  'plots': MultiselectInput,
  'color_order': dropdownInput,
  'order_when_continuous_color': checkboxInput,
  'size': sliderInput,
  'scale_by': dropdownInput
}