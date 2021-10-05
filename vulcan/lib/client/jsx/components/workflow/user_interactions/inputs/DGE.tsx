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
import { flattenStringOptions, joinNesting, StringOptions } from './monoids';
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
  const allData = useMemoized(joinNesting, data);
  const options = Object.keys(allData);
  const subOptions: DataEnvelope<StringOptions> = useMemo(() => {
    let opts: DataEnvelope<StringOptions> = {}
    Object.entries(allData).map( ([key,value]) => {
      const vals=Object.values(value as string);
      // String values only! 
      if (_.every(vals, _.isString)) {
        opts[key] = [...new Set(vals)];
      }
    })
    return opts;
  }, [data])
  
  const [doSubset, setSubset] = useState(false);
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
  function toggleSubsetting() {
    if (doSubset) removeSubset();
    if (!doSubset) addSubset();
    setSubset(!doSubset)
  }
  
  const setDEMethod = (method: string, vals = value) => {
    let newVals = output_sets[method]
    if (doSubset) {
      newVals['subset_meta'] = vals['subset_meta']
      newVals['subset_use'] = vals['subset_use']
    }
    onChange(some(newVals));
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
  
  console.log(value)
  
  return (
    <div>
      Step 1: Select your DE Question Type:
      <DropdownAutocomplete
          list={Object.keys(method_labels)}
          value={value['method']}
          onSelect={(val: string) => setDEMethod(method_labels[val])}
        />
      
      <hr/> 
      {"Use all cells or just a subset? "}
      <Button
          variant="outlined"
          color="primary"
          onClick={() => {toggleSubsetting()}}
          >
          Toggle Subsetting
        </Button>
      {SubsetComps(value, options, subOptions, updateValue)}
      
      <hr/>
      {DEComps(value, options, subOptions, updateValue)}
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
    'de_meta': null
  },
  'btwn-sets': {
    'method': 'btwn-sets',
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

const SubsetComps = (vals: DataEnvelope<any>, opts: string[], subopts: DataEnvelope<StringOptions>, changeFxn: Function) => {
  let value_select = null;
  if (Object.keys(vals).includes('subset_meta')) {
    if (vals['subset_meta']!=null) {
      value_select = MultiselectInput(
        'subset_use', changeFxn, vals['subset_use'],
        'Labels to keep', subopts[(vals['subset_meta'])] as string[])
    }
    return(
      <div>
        {dropdownInput(
          'subset_meta', changeFxn, vals['subset_meta'],
          "Subset by:", opts)}
        {value_select}
      </div>
    )
  }
  return(<div></div>)
}

const DEComps = (vals: DataEnvelope<any>, opts: string[], subopts: DataEnvelope<StringOptions>, changeFxn: Function) => {
  console.log(Object.keys(vals));
  let value_select_1 = null;
  let value_select_2 = null;
  if (Object.keys(vals).includes('de_meta')) {
    if (vals['de_meta']!=null) {
      if (Object.keys(vals).includes('de_group_1')) {
        value_select_1 = MultiselectInput(
          'de_group_1', changeFxn, vals['de_group_1'],
          'Labels to for Group-1', subopts[(vals['de_meta'])] as string[])
      }
      if (Object.keys(vals).includes('de_group_2')) {
        value_select_2 = MultiselectInput(
          'de_group_2', changeFxn, vals['de_group_2'],
          'Labels to keep', subopts[(vals['de_meta'])] as string[])
      }
    }
    return(
      <div>
        {dropdownInput(
          'de_meta', changeFxn, vals['de_meta'],
          'DiffEexp by', opts)}
        {value_select_1}
        {value_select_2}
      </div>
    )
  }
  return(<div></div>)
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