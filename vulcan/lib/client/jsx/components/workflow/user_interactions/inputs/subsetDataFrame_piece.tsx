import React, { useMemo } from 'react';
import {DataEnvelope} from './input_types';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import { dropdownInput, MultiselectInput, rangeInput, val_wrap } from './user_input_pieces';
import { flattenStringOptions } from './monoids';

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
      
      const updateColDef = (newColDef: (string|number|null)[], index: number, current_full = values) => {
        let next_full = {...current_full};
        next_full['methods'][index] = newColDef
        return next_full
      }
      
      const single_definition = (index: number, Def: (string|number|null)[] = [null]) => {
          
        const updateDef = (newDef: (string|number|null)[], x: string) => {
          // updateFxn for inner calls which will not include the column_name
          const col = [Def[0]]
          changeFxn(updateColDef(col.concat(...newDef), index), key)
        }
        
        const base = (
          <div key={key+index}>
            Target:
            <DropdownAutocomplete
              list={columns}
              value={Def[0]}
              onSelect={(newVal: string) => changeFxn(updateColDef([newVal], index), key)}
              sorted={sorted}
            />
          </div>)
        
        if (Def[0]==null) return base;
        const target_data = Object.values(data_frame[Def[0]])
        
        let inner_def = [...Def]
        inner_def.shift()
        
        let def_comps = null;
                      
        if (typeof target_data[0] == "number") {
          if (inner_def.length == 0) {
            inner_def = ["exactly", null, "below", null]
          }
          def_comps = rangeInput(
            key+index+index, updateDef, inner_def,
            "Values to target:"
            )
        }
        
        if (typeof target_data[0] == "string") {
          function onlyUnique(value: any, index: number, self: any) {
            return self.indexOf(value) === index;
          }
          const options = target_data.filter(onlyUnique) as string[];
          def_comps = MultiselectInput(
            key+index+index, updateDef, inner_def as string[],
            "Values to target:", options
            )
        }
        
        if (typeof target_data[0] == "boolean") {
          let use = ""
          if (inner_def.length > 0) {
            use = inner_def[0]
          }
          def_comps = dropdownInput(
            key+index+index, updateDef, use,
            "Values to target:", ["True", "False"], false
            )
        }
        
        return(
          <div key={key+index}>
            {base}
            {def_comps}
          </div>
        )
      }
      
      return(
        <div key={key}>
          {label}
          {values['methods'].map((Def, index) => {
            return single_definition(index, Def)
          })}
          Add another?
        </div>
      )
    }