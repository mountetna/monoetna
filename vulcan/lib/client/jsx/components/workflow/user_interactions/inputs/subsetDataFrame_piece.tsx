import React, { useMemo } from 'react';
import {DataEnvelope} from './input_types';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import { checkboxInput, dropdownInput, MultiselectInput, rangeInput } from './user_input_pieces';
import { Button } from '@material-ui/core';

/*
This script defines a component that behaves like all other 'user_input_pieces'.
It's a piece that recieves, as additional inputs:
 - 'data_frame' = a hash that corresponds to the json file of python pandas.DataFrame. Keys are column names, values are arrays of data where all elements have the same type.
 - 'sorted' = whether column name dropdowns should be alphanumerically sorted versus left in the same order as they are in the 'data_frame'.
*/


export function subsetDataFrameInput(
    key: string = "filler", changeFxn: Function, value: DataEnvelope<(string|number|null)[][]> = {'methods':[[null]]},
    label: string, data_frame: DataEnvelope<any>, sorted = false) {
      
      const values = {...value}
      const columns: string[] = Object.keys(data_frame)
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
      
      const single_definition = (index: number, Def: (string|number|null)[] = [null]) => {
        
        // updateFxns for inner piece calls
        // 'x' is a dummy variable needed because each piece will send a dummy key in addition the 'newDef' value.
        const updateLogic = (newLogic: string|null, x: string) => {
          updateCurrent(null, newLogic, index)
        }
        const updateDef = (newDef: (string|number|null)[], x: string) => {
          // Inners will not have the column_name
          const col = [Def[0]]
          updateCurrent(col.concat(...newDef), null, index)
        }
        
        const clearDef = () => {
          const full = {...values};
          // methods
          let next_methods = full['methods']
          next_methods.splice(index,1)
          if (next_methods.length==0) {next_methods = [[null]]}
          // logic
          changeFxn(
            {methods: next_methods},
            key)
        }
        
        let this_def_comps = null;
        let def_comp = null;
        const pick_column = (
          <div key={key+index} style={{display: 'inline-flex'}}>
            Feature {" "+(index+1)}:
            <DropdownAutocomplete
              list={columns}
              value={Def[0]}
              onSelect={(newVal: string) => updateCurrent([newVal], null, index)}
              sorted={sorted}
            />
          </div>
        )
        const clear_comp = (
          <Button onClick={(x) => {clearDef()}}>
            Clear
          </Button>
        )
        const logic_comp = (index == 0) ? null : 
          dropdownInput(
              key+index+index+index, updateLogic, values['logic'][index-1][0] as string,
              "Combination Logic:", ["and", "or"], false
            )
        
        if (Def[0]!=null) {

          const target_data = Object.values(data_frame[Def[0]])
          
          let inner_def = [...Def]
          inner_def.shift()
                        
          if (typeof target_data[0] == "number") {
            if (inner_def.length == 0) {
              inner_def = ["exactly", null, "below", null]
            }
            def_comp = rangeInput(
              key+index+index, updateDef, inner_def,
              ""
              )
          }
          
          if (typeof target_data[0] == "string") {
            function onlyUnique(value: any, index: number, self: any) {
              return self.indexOf(value) === index;
            }
            const options = target_data.filter(onlyUnique) as string[];
            def_comp = MultiselectInput(
              key+index+index, updateDef, inner_def as string[],
              "Keep:", options
              )
          }
          
          if (typeof target_data[0] == "boolean") {
            let use = ""
            if (inner_def.length > 0) {
              use = inner_def[0] as string
            }
            def_comp = dropdownInput(
              key+index+index, updateDef, use,
              "", ["True", "False"], false
              )
          }
        }
        
        this_def_comps = (
          <div style={{display: 'inline-flex'}}>
            <div>
              {logic_comp}
              {pick_column}
              {def_comp}
            </div>
            {clear_comp}
          </div>
        )
        
        return(
          <div key={key+index} style={{
            paddingLeft: "5px",
            paddingTop: "10px"}}>
            {this_def_comps}
          </div>
        )
      }
      
      const startOrClear = (doSubset: boolean, x:any) => {
        let new_full = {};
        if (doSubset) {
          new_full = {
            'methods':[[null]],
            'logic':[[]]
          }
        }
        changeFxn(new_full, key)
      }
      
      const base = checkboxInput(
        key, startOrClear, Object.keys(values).includes('methods'),
        label)
      
      let inner = null;
      if (values && values['methods']) {
        inner = (
          <div>
            {values['methods'].map((Def, index) => {
              return single_definition(index, Def)
              })
            }
            <Button onClick={(x) => {updateCurrent([null], "and", values['methods'].length)}}>
              Add another condition?
            </Button>
          </div>
        )
      }
      return(
        <div key={key}>
          <p></p>
          {base}
          {inner}
          <p></p>
        </div>
      )
    }
