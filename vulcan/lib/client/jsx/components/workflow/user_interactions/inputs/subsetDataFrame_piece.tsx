import React, { useMemo } from 'react';
import {DataEnvelope} from './input_types';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import { dropdownInput, MultiselectInput, rangeInput } from './user_input_pieces';

/*
This script defines a component that behaves like all other 'user_input_pieces'.
It's a piece that recieves, as additional inputs:
 - 'data_frame' = a hash that corresponds to the json file of python pandas.DataFrame. Keys are column names, values are arrays of data where all elements have the same type.
 - 'sorted' = whether column name dropdowns should be alphanumerically sorted versus left in the same order as they are in the 'data_frame'.
*/


export function subsetDataFrameInput(
    key: string = "filler", changeFxn: Function, value: DataEnvelope<(string|number|null)[][]> = {'methods':[[null]]},
    label: string, data_frame: DataEnvelope<any>, sorted = false) {
      
      const columns: string[] = Object.keys(data_frame)
      
      const updateColDef = (newColDef: (string|number|null)[], index: number, current_full = value) => {
        let next_full = {...current_full};
        next_full['methods'][index] = newColDef
        return next_full
      }
      
      const single_definition = (index: number, Def: (string|number|null)[] = [null]) => {
          
        const updateDef = (newDef: (string|number|null)[], key: string) => {
          // updateFxn for inner calls which will not include the column_name
          const col = Def[0]
          updateColDef([col].concat(...newDef), index)
        }
        
        const Definition = (target: string|null) => {
          if (target==null) return null
          const target_data = data_frame[target]
            
          if (typeof target_data[0] == "number") {
            return rangeInput(
              key+index+index, updateDef, ["exactly", null, "below", null],
              "Values to target:"
              )
          }
          
          if (typeof target_data[0] == "string") {
            return MultiselectInput(
              key+index+index, updateDef, null,
              "Values to target:", columns
              )
          }
          
          if (typeof target_data[0] == "boolean") {
            return dropdownInput(
              key+index+index, updateDef, null,
              "Values to target:", ["True", "False"], false
              )
          }
          
          return (
              <div> Functionality Under Development </div>
          )
        }
        
        return(
          <div key={key+index}>
            Target:
            <DropdownAutocomplete
              list={columns}
              value={Def[0]}
              onSelect={(newVal: string) => changeFxn(updateColDef([newVal], index), key)}
              sorted={sorted}
            />
            {Definition(Def[0] as string|null)}
          </div>
        )
      }
      
      return(
        <div key={key}>
          {label}
          {value['methods'].map((Def, index) => {
            return single_definition(index, Def)
          })}
          Add another?
        </div>
      )
    }