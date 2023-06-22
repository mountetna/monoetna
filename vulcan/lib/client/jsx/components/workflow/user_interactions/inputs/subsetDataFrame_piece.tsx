import React, {useMemo} from 'react';
import {DataEnvelope} from './input_types';
import {
  checkboxPiece,
  dropdownPiece,
  key_wrap,
  MultiselectPiece,
  nestedDropdownPiece,
  rangePiece,
  arrayLevels
} from './user_input_pieces';
import {Button, PropTypes} from '@material-ui/core';
import DeleteIcon from '@material-ui/icons/Delete';
import { nestedOptionSet } from './visualizations';

/*
This script defines a component that behaves like all other 'user_input_pieces'.
See 'user_input_pieces.tsx' for additional documentation.
It's a piece that recieves, as additional inputs:
 - 'data_summary' = a hash that corresponds to the json file of python pandas.DataFrame, OR a summary of unique values for all discrete or boolean data from a data_frame or similar data structure. Keys are column names, values are arrays of data where all elements have the same type.
 - 'additional_numeric_options' = a string array naming any numeric data columns that can be subset on but are not included in the 'data_summary'.
 - 'sorted' = whether column name dropdowns should be alphanumerically sorted versus left in the same order as they are in the 'data_frame'.
Conceptually, it's easiest to think about this piece as if data_summary represents a data_frame and this component is being used to define any number of filters that each target a subsets of values a single column of that data_frame.
For practical purposes (specifically: 20k+ genes in genomic datasets (-> data_frame columns), or 100k+ cells (-> data_frame rows)),
the component was abstracted so that 'data_summary' could:
 - exclude numeric 'columns'
 - recieve only the unique values of discrete 'columns', rather than requiring ALL data points.
*/

const emptyMethod = [null];

function targetSelectionComponent(
  def: (string | number | null)[],
  data_summary: DataEnvelope<any>,
  updateFxn: Function,
  key: string,
  index: number
) {
  if (def[0] == emptyMethod[0]) return null;

  let inner_def = [...def];
  inner_def.shift();

  if ( ! Object.keys(data_summary).includes(def[0] as string)) {
    // Case = subsetting on data from 'additional_numeric_options'
    return rangePiece(
      key + index + index,
      updateFxn,
      inner_def.length > 0 ? inner_def : undefined,
      'Values to keep'
    );
  }

  const target_data = Object.values(data_summary[def[0] as string]);

  if (typeof target_data[0] == 'number') {
    return rangePiece(
      key + index + index,
      updateFxn,
      inner_def.length > 0 ? inner_def : undefined,
      'Values to keep'
    );
  }
  if (typeof target_data[0] == 'string') {
    // Need to build the unique options set first.
    const options = arrayLevels(target_data) as string[];
    return MultiselectPiece(
      key + index + index,
      updateFxn,
      inner_def as string[],
      'Values to keep',
      options
    );
  }
  if (typeof target_data[0] == 'boolean') {
    return dropdownPiece(
      key + index + index,
      updateFxn,
      inner_def.length > 0 ? (inner_def[0] as string) : undefined,
      'Value to keep',
      ['True', 'False'],
      false
    );
  }
  return (
    <div>
      Not yet implemented Data Taype. Please followup with someone from the Data
      Library Team!
    </div>
  );
}

const singleMethod = (
  def: (string | number | null)[] = emptyMethod,
  index: number,
  data_summary: DataEnvelope<any>,
  subset_options: DataEnvelope<DataEnvelope<DataEnvelope<null>|null>|null> | string[],
  updateCurrent: Function,
  overallChangeFxn: Function,
  key: string,
  values: DataEnvelope<(string | number | null)[][]>,
  sorted: boolean,
  color: PropTypes.Color
) => {
  // updateFxns for inner piece calls
  // 'x' is a dummy variable needed because each piece will send a dummy key in addition the 'newDef' value.
  const updateLogic = (newLogic: string | null, x: string) => {
    updateCurrent(null, newLogic, index);
  };
  const updateDefTargets = (newDef: (string | number | null)[], x: string) => {
    // Data target inners will not have the column_name
    const col = [def[0]];
    updateCurrent(col.concat(...newDef), null, index);
  };
  const updateDefColumn = (newCol: string, x: string) => {
    // (Resets away any data targets)
    updateCurrent([newCol], null, index);
  };
  
  const clearDef = () => {
    const full = {...values};
    // methods (Current at index, must remain/become [emptyMethod] or lose current)
    let next_methods = full['methods'];
    next_methods.splice(index, 1);
    if (next_methods.length == 0) {
      next_methods = [emptyMethod];
    }
    // logic (Current at index-1, must remain/become [[]], lose current, or lose next if clicked for 1st def)
    let next_logic = full['logic'];
    if (next_logic.length == 1) {
      next_logic = [[]];
    } else {
      if (index == 0) {
        next_logic.splice(index, 1);
      }
      if (index > 0) {
        next_logic.splice(index - 1, 1);
      }
    }
    overallChangeFxn({methods: next_methods, logic: next_logic}, key);
  };

  const pick_column = Array.isArray(subset_options) ? dropdownPiece(
    key + index,
    updateDefColumn,
    def[0] as string,
    'Condition ' + (index + 1) + ', Feature',
    subset_options,
    sorted
  ) : nestedDropdownPiece(
    key + index,
    updateDefColumn,
    def[0] as string,
    'Condition ' + (index + 1) + ', Feature',
    subset_options
  );
  const clear_comp = (
    <Button
      color={color}
      onClick={(x) => {
        clearDef();
      }}
    >
      <DeleteIcon fontSize='small' />
    </Button>
  );
  const logic_comp =
    index == 0
      ? null
      : dropdownPiece(
          key + index + index + index,
          updateLogic,
          values['logic'][index - 1][0] as string,
          'Combination Logic',
          ['and', 'or'],
          false
        );
  const def_comp = targetSelectionComponent(
    def,
    data_summary,
    updateDefTargets,
    key,
    index
  );

  return (
    <div
      key={key + index}
      style={{
        display: 'inline-flex'
      }}
    >
      <div>
        {logic_comp}
        {pick_column}
        <div style={{paddingLeft: 10}}>{def_comp}</div>
      </div>
      {clear_comp}
    </div>
  );
};

export function subsetDataFramePiece(
  key: string = 'filler',
  changeFxn: Function,
  value: DataEnvelope<(string | number | null)[][]> = {},
  label: string,
  data_summary_ori: DataEnvelope<any>,
  sorted = false,
  color: PropTypes.Color = 'primary',
  additional_numeric_options: nestedOptionSet | null = null,
  data_summary_options_label: string = 'data_frame columns'
) {
  const values = {...value};
  const data_summary = {...data_summary_ori};
  // Put together all subsetting options
  let subset_options: string[] | nestedOptionSet = Object.keys(data_summary);
  if (additional_numeric_options!=null) {
    subset_options = {
      ...additional_numeric_options,
      [data_summary_options_label]: key_wrap([...subset_options])}
  }

  const startOrClear = (doSubset: boolean, x: any) => {
    const new_full = doSubset
      ? {
          methods: [emptyMethod],
          logic: [[]]
        }
      : {};
    changeFxn(new_full, key);
  };
  const updateCurrent = (
    newColDef: (string | number | null)[] | null,
    newLogic: string | null,
    index: number,
    current_full = values
  ) => {
    let next_full = {...current_full};
    if (newColDef) {
      next_full['methods'][index] = newColDef;
    }
    if (newLogic && index > 0) {
      next_full['logic'][index - 1] = [newLogic];
    }
    changeFxn(next_full, key);
  };

  const base = checkboxPiece(
    key,
    startOrClear,
    Object.keys(values).includes('methods'),
    label
  );

  const meat =
    values && values['methods'] ? (
      <div
        style={{
          paddingLeft: '15px',
          paddingTop: '2px'
        }}
      >
        {values['methods'].map((def, index) => {
          return singleMethod(
            def,
            index,
            data_summary,
            subset_options,
            updateCurrent,
            changeFxn,
            key,
            values,
            sorted,
            color
          );
        })}
        <br></br>
        <Button
          color={color}
          onClick={(x) => {
            updateCurrent(emptyMethod, 'and', values['methods'].length);
          }}
        >
          Add another condition?
        </Button>
      </div>
    ) : null;
  // console.log(values);
  return (
    <div key={key}>
      {base}
      {meat}
    </div>
  );
}
