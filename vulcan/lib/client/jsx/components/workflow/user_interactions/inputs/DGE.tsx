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
import { flattenStringOptions, joinNesting } from './monoids';
import { useMemoized } from '../../../../selectors/workflow_selectors';

/*
This input is closely tied to a python script for running scanpy-based differential expression.

Input:
- The metadata data.frame, as json, of an single-cell object.

Output:
- A dictionary containing:
    - method = a string of "btwn-all-de-groups", "btwn-sets", "btwn-sets-multiple-groups"
    - (optional) subset_meta = a string naming a metadata to subset by (when: subsetting needed)
    - (optional) subset_use = a string vector of values of the subset meta to keep (when: subsetting needed)
    - de_meta = a string naming a metadata to compare between
    - de_group_1 = a string vector of values of the de_meta to either compare among (when: "btwn-all-de-groups") or use as group1 (when:"btwn-sets", "btwn-sets-multiple-groups")
    - (optional) de_group_2 = a string vector of values of the de_meta to use as group2 (when: "btwn-sets", "btwn-sets-multiple-groups")
    - (optional) groups_meta = a string naming a metadata to use as cell groupings (when: "btwn-sets-multiple-groups")
    - (optional) groups_use = a string vector indicating which values of the groupings meta to investigate within (when: "btwn-sets-multiple-groups")

Labels:
- Between all labels of a particular metadata
- Between particular labels
- Between particular labels in multiple cell groups
*/

export default function DiffExpSC({
  data, onChange, ...props
}: WithInputParams<{}, DataEnvelope<any>, any>, plotType: string) {
  const value = useSetsDefault({method: null}, props.value, onChange);
  const allOptions = useMemoized(joinNesting, data);
  const options = Object.keys(allOptions);
  const [doSubset, setSubset] = useState(false);
  
  /*
  const subOptions: DataEnvelope<string[]> = useMemo(() => {
    if (data == null) return [];
    const vals = {...data}
    let opts = {...data}
    Object.entries(vals).map( ([key,value]) => {
      opts[key] = flattenStringOptions({'0': value})
    })
    return opts;
  }, [data])
  */
  const addSubset = (vals = value) => {
    vals['subset_meta'] = null
    vals['subset_use'] = null
    onChange(some(vals));
  }
  
  const removeSubset = (vals = value) => {
    if (Object.keys(vals).includes('subset_meta')) delete vals['subset_meta']
    if (Object.keys(vals).includes('subset_use')) delete vals['subset_use']
    onChange(some(vals));
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
 
  const SubsetComps = (vals = value, opts = options, allopts = allOptions) => {
    console.log(Object.keys(vals));
    let value_select = null;
    if (Object.keys(vals).includes('subset_meta')) {
      if (vals['subset_meta']!=null) {
        value_select = MultiselectInput(
          'subset_use', updateValue, vals['subset_use'],
          'Labels to keep', allopts[vals['subset_meta']] as string[])
      }
      return(
        <div>
          {dropdownInput(
            'subset_meta', updateValue, vals['subset_meta'],
            "Subset by:", opts)}
          {value_select}
        </div>
      )
    } else {
      return(<div></div>)
    }
  }
  
  function toggleSubsetting(vals = value) {
    if (doSubset) removeSubset();
    if (!doSubset) addSubset();
    setSubset(!doSubset)
  }
  
  console.log(options)
  console.log(value)
  
  return (
    <div>
      Differential Expression Question:
      <DropdownAutocomplete
          list={Object.keys(method_labels)}
          value={value['method']}
          onSelect={(val: string) => updateValue(method_labels[val], 'method')}
        />
      Use all cells or just a subset?
      <Button
          variant="contained"
          color="primary"
          onClick={() => {toggleSubsetting()}}
          >
          Subsetting
        </Button>
      {SubsetComps()}
      {/* {DEComps} */}
      {/* {GroupsComps} */}
    </div>
  );

};

const method_labels: DataEnvelope<string> = {
  '1 btwn-all-de-groups: Between All labels of one metadata': 'btwn-all-de-groups',
  '2 btwn-sets: Between particular labels': 'btwn-sets',
  '3 btwn-sets-multiple-groups: Between particular labels, in multiple cell groups': 'btwn-sets-multiple-groups'
}

const output_sets: DataEnvelope<DataEnvelope<string|string[]|null>> = {
  'btwn-all-de-groups': {
    'method': 'btwn-all-de-groups',
    'de_meta': null,
    'de_group_1': null
  },
  'btwn-sets': {
    'method': 'btwn-all-de-groups',
    'de_meta': null,
    'de_group_1': null,
    'de_group_2': null
  },
  'btwn-sets-multiple-groups': {
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