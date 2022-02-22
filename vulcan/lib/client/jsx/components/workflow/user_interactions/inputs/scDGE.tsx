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
import { val_wrap, MultiselectPiece, dropdownPiece } from './user_input_pieces';
import { subsetDataFramePiece } from './subsetDataFrame_piece';
import { Button } from '@material-ui/core';

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
  const value = useSetsDefault({method: null, subset: {}}, props.value, onChange);
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
  
  const setDEMethod = (label: string, x: string, label_tags = method_labels, vals = {...value}) => {
    // x = filler for 'key' provided by input piece which will not match the current aims.
    let newVals = {...output_sets[label_tags[label]]}
    newVals['subset'] = vals['subset']
    onChange(some(newVals));
  }
  
  const updateValue = (newValue: any, key: string, prevValues = {...value}) => {
    prevValues[key] = newValue;
    onChange(some(prevValues));
  };

  const warn =
    <div>
      <p></p>
      IMPORTANT NOTE: The scanpy DE algorithm cannot correct for batch effects. We are currently investigating alternatives to add such functionality.
      <p></p>
    </div>

  // Bringing pre-subsetDataFrameInput value-structure up-to-date.
  let SubsetComps = null;
  if (value) {
    if (value['subset']) {
    SubsetComps = subsetDataFramePiece(
      "subset", updateValue, value['subset'], "Step 2: Subset by features to focus on certain cells? ",
      {...allData}, false, "secondary")
    } else {
      const updateOldSubset = () => {
        let new_subset = ( !(value['subset_meta'] && value['subset_use'])) ? {} : {
          'methods': [ [value['subset_meta']].concat(value['subset_use']) ],
          'logic': [[]]
        }
        let prevValues = {...value}
        delete prevValues['subset_meta']
        delete prevValues['subset_use']
        updateValue(new_subset, 'subset', prevValues)
      }
      SubsetComps = 
        <Button
          color={"primary"}
          onClick={(x) => {updateOldSubset()}}>
          Step 2: Subset... UPDATE NEEDED. Click to update.
        </Button>
    }
  }

  return (
    <div>
      {dropdownPiece(
        "method", setDEMethod, value['method'],
        "Step 1: Select your DE Question Type: ",
        Object.keys(method_labels), true)}
      <hr/>
      {SubsetComps}
      {DEComps(value, options, updateValue)}
      {GroupComps(value, options, updateValue)}
      {warn}
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
          value_select_1 = MultiselectPiece(
            'de_group_1', changeFxn, vals['de_group_1'],
            'Labels to include in Set-1', opts[(vals['de_meta'])] as string[])
        }
        if (Object.keys(vals).includes('de_group_2')) {
          value_select_2 = MultiselectPiece(
            'de_group_2', changeFxn, vals['de_group_2'],
            'Labels to include in Set-2', opts[(vals['de_meta'])] as string[])
        }
      }
      
      return(
        <div>
          <hr/>
          {"Step 3: What labels do you want to compare? "}
          {dropdownPiece(
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
          value_select = MultiselectPiece(
            'group_use', changeFxn, vals['group_use'],
            'Labels to target', opts[(vals['group_meta'])] as string[])
        }
      }
      
      return(
        <div>
          <hr/>
          {"Step 4: Within what comparison groups?"}
          {dropdownPiece(
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
