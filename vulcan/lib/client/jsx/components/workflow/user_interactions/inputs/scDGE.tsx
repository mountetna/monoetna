// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useMemo, useState} from 'react';
import * as _ from 'lodash';

import {DataEnvelope, WithInputParams} from './input_types';
import { useSetsDefault } from './useSetsDefault';
import { maybeOfNullable, some, withDefault, Maybe } from '../../../../selectors/maybe';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import BooleanInput from './boolean';
import { joinNesting, StringOptions } from './monoids';
import { useMemoized } from '../../../../selectors/workflow_selectors';
import { val_wrap, MultiselectInput, dropdownInput } from './user_input_pieces';

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
    - (optional) de_group_1 = a string vector of values of the de_meta to use as group1 (when: "btwn-sets", "btwn-sets-multiple-groups")
    - (optional) de_group_2 = a string vector of values of the de_meta to use as group2 (when: "btwn-sets", "btwn-sets-multiple-groups")
    - (optional) groups_meta = a string naming a metadata to use as cell groupings (when: "btwn-sets-multiple-groups")
    - (optional) groups_use = a string vector indicating which values of the groupings meta to investigate within (when: "btwn-sets-multiple-groups")
*/

export default function DiffExpSC({
  data, onChange, ...props
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  const value = useSetsDefault({method: null}, props.value, onChange);
  const allData = useMemoized(joinNesting, data);
  const options: DataEnvelope<StringOptions> = useMemo(() => {
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
  const addSubset = (vals = {...value}) => {
    vals['subset_meta'] = null
    vals['subset_use'] = null
    onChange(some(vals));
  }
  const removeSubset = (vals = {...value}) => {
    if (Object.keys(vals).includes('subset_meta')) delete vals['subset_meta']
    if (Object.keys(vals).includes('subset_use')) delete vals['subset_use']
    onChange(some(vals));
  }
  function toggleSubsetting(subset: boolean) {
    if (!subset) removeSubset();
    if (subset) addSubset();
    setSubset(subset)
  }
  
  const setDEMethod = (method: string, vals = {...value}) => {
    let newVals = {...output_sets[method]}
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
  
  return (
    <div>
      Step 1: Select your DE Question Type:
      <DropdownAutocomplete
          list={Object.keys(method_labels)}
          value={value['method']}
          onSelect={(val: string) => setDEMethod(method_labels[val])}
        />
      
      {SubsetComps(value, options, updateValue, doSubset, toggleSubsetting)}
      {DEComps(value, options, updateValue)}
      {GroupComps(value, options, updateValue)}
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
    'method': 'btwn-sets-multiple-groups',
    'de_meta': null,
    'de_group_1': null,
    'de_group_2': null,
    'group_meta': null,
    'group_use': null
  }
}

const SubsetComps = (vals: DataEnvelope<any>, opts: DataEnvelope<StringOptions>, changeFxn: Function, doSubset: boolean, toggleSubsetting: any) => {
  
  const comps = useMemo(() => {
    const base =
      <div>
        <hr/> 
        {"Step 2: Focus on / remove certain cells? "}
        <div>
          <BooleanInput
            key='doSubset'
            label='Subset'
            value={maybeOfNullable(doSubset)}
            data={val_wrap(doSubset)}
            onChange={ value => toggleSubsetting(withDefault(value, false))}
          />
        </div>
      </div>
    
    let value_select = null;
    if (Object.keys(vals).includes('subset_meta')) {
      if (Object.keys(opts).length>0 && vals['subset_meta']!=null) {
        value_select = MultiselectInput(
          'subset_use', changeFxn, vals['subset_use'],
          'Labels to keep', opts[(vals['subset_meta'])] as string[])
      }
      return(
        <div>
          {base}
          {dropdownInput(
            'subset_meta', changeFxn, vals['subset_meta'],
            "Subset by:", Object.keys(opts))}
          {value_select}
        </div>
      )
    }
    return(base)
  }, [vals, opts])
  
  return comps
}

const DEComps = (
  vals: DataEnvelope<any>,
  opts: DataEnvelope<StringOptions>,
  changeFxn: Function) => {
  
  const comps = useMemo(() => {
    if (Object.keys(vals).includes('de_meta')) {
      
      let value_select_1 = null;
      let value_select_2 = null;
      if (Object.keys(opts).length>0 && vals['de_meta']!=null) {
        if (Object.keys(vals).includes('de_group_1')) {
          value_select_1 = MultiselectInput(
            'de_group_1', changeFxn, vals['de_group_1'],
            'Labels to for Group-1', opts[(vals['de_meta'])] as string[])
        }
        if (Object.keys(vals).includes('de_group_2')) {
          value_select_2 = MultiselectInput(
            'de_group_2', changeFxn, vals['de_group_2'],
            'Labels to keep', opts[(vals['de_meta'])] as string[])
        }
      }
      
      return(
        <div>
          <hr/>
          {"Step 3: What labels do you want to compare? "}
          {dropdownInput(
            'de_meta', changeFxn, vals['de_meta'],
            'DiffExp by:', Object.keys(opts))}
          {value_select_1}
          {value_select_2}
        </div>
      )
    }
    return(<div></div>)
  }, [vals, opts])
  
  return comps
}

const GroupComps = (
  vals: DataEnvelope<any>,
  opts: DataEnvelope<StringOptions>,
  changeFxn: Function) => {

  const comps = useMemo(() => {
    if (Object.keys(vals).includes('group_meta')) {
      
      let value_select = null;
      if (Object.keys(opts).length>0 && vals['group_meta']!=null) {
        if (Object.keys(vals).includes('group_use')) {
          value_select = MultiselectInput(
            'group_use', changeFxn, vals['group_use'],
            'Labels to focus on', opts[(vals['group_meta'])] as string[])
        }
      }
      
      return(
        <div>
          <hr/>
          {"Step 4: Between what comparison groups?"}
          {dropdownInput(
            'group_meta', changeFxn, vals['group_meta'],
            'Group by:', Object.keys(opts))}
          {value_select}
        </div>
      )
    }
    return(<div></div>)
  }, [vals, opts])
  
  return comps
}
