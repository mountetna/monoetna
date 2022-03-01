import React, { useMemo } from 'react';
import {DataEnvelope} from './input_types';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import { checkboxPiece, dropdownPiece, MultiselectPiece, rangePiece } from './user_input_pieces';
import { Button, PropTypes } from '@material-ui/core';
import DeleteIcon from '@material-ui/icons/Delete';

/*
This script defines a component that behaves like all other 'user_input_pieces'.
See 'user_input_pieces.tsx' for additional documentation.
It's a piece that recieves, as additional inputs:
 - 'data_frame' = a hash that corresponds to the json file of python pandas.DataFrame. Keys are column names, values are arrays of data where all elements have the same type.
 - 'sorted' = whether column name dropdowns should be alphanumerically sorted versus left in the same order as they are in the 'data_frame'.
*/

export function subsetDataFramePiece(
  key: string = "filler", changeFxn: Function, value: DataEnvelope<(string|number|null)[][]> = {},
  label: string, data_frame_ori: DataEnvelope<any>, sorted = false, color: PropTypes.Color = "primary") {
    
  const values = {...value}
  const data_frame = {...data_frame_ori}
  
  const startOrClear = (doSubset: boolean, x:any) => {
    const new_full = (doSubset) ? {
      'methods':[emptyMethod],
      'logic':[[]]
    } : {}
    changeFxn(new_full, key)
  }
  const updateCurrent = (
    newColDef: (string|number|null)[] | null,
    newLogic: string | null,
    index: number,
    current_full = values
    ) => {
    let next_full = {...current_full};
    if (newColDef) {
      next_full['methods'][index] = newColDef
    }
    if (newLogic && index>0) {
      next_full['logic'][index-1] = [newLogic]
    }
    changeFxn(next_full, key)
  }
  
  const base = checkboxPiece(
    key, startOrClear, Object.keys(values).includes('methods'),
    label)
  
  const meat = (values && values['methods']) ? (
    <div style={{
      paddingLeft: "5px",
      paddingTop: "2px",}}>
      {values['methods'].map((def, index) => {
        return singleMethod(
          def, index, data_frame, updateCurrent, changeFxn, key, values, sorted, color)
        })
      }
      <br></br>
      <Button
        color={color}
        onClick={(x) => {updateCurrent(emptyMethod, "and", values['methods'].length)}}>
        Add another condition?
      </Button>
    </div>
    ) : null;
  // console.log(values);
  return(
    <div key={key}>
      {base}
      {meat}
    </div>
  )
}

const emptyMethod = [null]

const singleMethod = (
  def: (string|number|null)[] = emptyMethod,
  index: number,
  data_frame: DataEnvelope<any>,
  updateCurrent: Function,
  overallChangeFxn: Function,
  key: string,
  values: DataEnvelope<(string|number|null)[][]>,
  sorted: boolean,
  color: PropTypes.Color
  ) => {
    
  // updateFxns for inner piece calls
  // 'x' is a dummy variable needed because each piece will send a dummy key in addition the 'newDef' value.
  const updateLogic = (newLogic: string|null, x: string) => {
    updateCurrent(null, newLogic, index)
  }
  const updateDefTargets = (newDef: (string|number|null)[], x: string) => {
    // Data target inners will not have the column_name
    const col = [def[0]]
    updateCurrent(col.concat(...newDef), null, index)
  }
  const updateDefColumn = (newCol: string, x: string) => {
    // (Resets away any data targets)
    updateCurrent([newCol], null, index)
  }
  
  const columns: string[] = Object.keys(data_frame)
  
  const clearDef = () => {
    const full = {...values};
    // methods (Current at index, must remain/become [emptyMethod] or lose current)
    let next_methods = full['methods']
    next_methods.splice(index,1)
    if (next_methods.length==0) {
      next_methods = [emptyMethod]
    }
    // logic (Current at index-1, must remain/become [[]], lose current, or lose next if clicked for 1st def)
    let next_logic = full['logic']
    if (next_logic.length == 1) {
      next_logic = [[]]
    } else {
      if (index==0) {
        next_logic.splice(index,1)
      }
      if (index>0) {
        next_logic.splice(index-1,1)
      }
    }
    overallChangeFxn(
      {methods: next_methods, logic: next_logic},
      key)
  }
  
  const pick_column = dropdownPiece(
    key+index, updateDefColumn, def[0] as string,
    "Condition " + (index+1) + ", Feature:", columns, sorted
  )
  const clear_comp = (
    <Button
      color={color}
      onClick={(x) => {clearDef()}}>
      <DeleteIcon fontSize="small"/>
    </Button>
  )
  const logic_comp = (index == 0) ? null : 
    dropdownPiece(
        key+index+index+index, updateLogic, values['logic'][index-1][0] as string,
        "Combination Logic:", ["and", "or"], false
      )
  const def_comp = targetSelectionComponent(
    def, data_frame, updateDefTargets, key, index)
  
  return(
    <div key={key+index}
      style={{
        display: 'inline-flex'}}>
      <div>
        {logic_comp}
        {pick_column}
        {def_comp}
      </div>
      {clear_comp}
    </div>
  )
}

function targetSelectionComponent(
  def: (string|number|null)[],
  data_frame: DataEnvelope<any>,
  updateFxn: Function,
  key: string,
  index: number
  ) {
    
  if (def[0]==emptyMethod[0]) return null
  
  let inner_def = [...def]
  inner_def.shift()
  
  const target_data = Object.values(data_frame[def[0] as string])
  
  if (typeof target_data[0] == "number") {
    return rangePiece(
      key+index+index, updateFxn,
      (inner_def.length > 0) ? inner_def : ["exactly", null, "below", null],
      ""
      )
  }
  if (typeof target_data[0] == "string") {
    // Need to build the unique options set first.
    function onlyUnique(value: any, index: number, self: any) {
      return self.indexOf(value) === index;
    }
    const options = target_data.filter(onlyUnique) as string[];
    return MultiselectPiece(
      key+index+index, updateFxn,
      inner_def as string[],
      "Keep:", options
      )
  }
  if (typeof target_data[0] == "boolean") {
    return dropdownPiece(
      key+index+index, updateFxn,
      (inner_def.length > 0) ? inner_def[0] as string : "",
      "Keep:", ["True", "False"], false
      )
  }
  return(
    <div>Not yet implemented Data Taype. Please followup with someone from the Data Library Team!</div>
  ) 
}
