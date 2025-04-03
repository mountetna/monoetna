
import React from 'react';
import {DataEnvelope} from '../../input_types';
import {
  dropdownPiece,
  MultiselectPiece,
  nestedDropdownPiece,
  rangePiece,
  arrayLevels,
  NumericConstraint,
  emptyNumericConstraintDef
} from './user_input_pieces';
import {Button, PropTypes} from '@material-ui/core';
import DeleteIcon from '@material-ui/icons/Delete';
import InputLabel from '@material-ui/core/InputLabel';
import { nestedOptionSet } from './utils';

/*
This script establishes various Pieces all aimed at data slicing for the purpose of extablishing data groupings.
Often, multiple groupings may be desired, so there are various flavors here.

For selection of a single group:
Often, selection of a grouping will rely on 1) a choice of data to select on, then 2) parameterization of selection contraints on that data.
Sometimes, selection of a grouping may require multiple such setups combined together with AND/OR logic.
E.g.:
  - 1) column "treatment" has value of "CD3-stim" or "CD3-CD28-stim", AND 2) column CD45 has value greater than 0. This might select cells with CD45-capture of stimulated samples. 

An individual constraint parameterization output style depends on the type of data:
  numeric:
    - length-4 array where elements are
      0: "exactly" OR "above"
      1: number OR null (where null means not yet selected)
      2: "exactly" OR "below"
      3: number OR null (where null means not yet selected)
    - Ex: ['exactly', 0, 'below', 5] ==> observations where the values is in [0,5) (or values >= 0 and < 5) are selected.
  string / categories:
    - any-length array of selected values
    - Ex: ["a", "c"] ==> If column contains ["a", "a", "b", "b", "c", "c"], the first and last 2 are selected.
  boolean:
    - length-1 array
    - Ex: ["true"] ==> observations where this data is True are selected.

When combining multiple data constraints into a single selection-definition, all constraints after the first will require a 'logic' element of 'AND' or 'OR'.

The overall format

In some cases, there may be desire to combine multiple selection methods.  E.g. 'condition' in ["A", "B"] & 'timepoint' in ["day0"].
This is not currently enabled here.

These pieces all assume there is a data_frame of data where:
- a 'data_summary' object ammounts to either the literal data_frame or a representation of it
  - in representation form: {[k: string]: string[]} where keys are column names of the dataframe that contain categorical data and values are either the literal values of the column, or the unique values of the original data.
- Assumption: all data of a given data column are of the same type
*/

/*
Data Provision:
Continuing with the way other vulcan UIs allow use of a smaller description of data-frame column contents rather than requiring the full data_frame as input:
- data_summary: DataEnvelope<any>, = the data_frame itself OR the set of unique values as string[] for any non-numeric columns of the actual data_frame.
- all_column_options?: nestedOptionSet | string[], = optional method for providing selection of data columns.
    - when data_summary receives the full data_frame itself, *not needed at all*, but can still be used to give some organization via catagories to column selection
    - when data_summary is just the summary of unique values of non-numeric columns, *absolutely needed*, provides names of all numeric columns that don't exist in the 'data_summary'.
*/

const emptyConstraintDef = [] as string[] | NumericConstraint;
type SingleSelectionConstraintDef = typeof emptyConstraintDef;
const emptySingleSelectionDefinition = {
  col: null as string | null,
  def: [] as SingleSelectionConstraintDef,
  logic: null as null | 'AND' | 'OR'
};
type SingleSelectionDefinition = typeof emptySingleSelectionDefinition;
export type SelectionDefinition = SingleSelectionDefinition[]
export const emptySelectionDefinition: SelectionDefinition = [emptySingleSelectionDefinition]

export function SingleConstraintDefPiece(
  key: string,
  changeFxn: (v: SingleSelectionConstraintDef, k?: string) => void,
  value: SingleSelectionConstraintDef,
  label: string = 'Targets',
  data_values: undefined | (string | number | boolean)[]
) {
  
  if ( !data_values || typeof data_values[0] == 'number') {
    // Special Case = subsetting on numeric data not described in the 'data_summary'
    return rangePiece(
      key,
      changeFxn,
      value.length > 0 ? value as NumericConstraint: emptyNumericConstraintDef,
      label
    );
  }
  if (typeof data_values[0] == 'string') {
    // Need to build the unique options set first.
    const options = arrayLevels(data_values) as string[];
    return MultiselectPiece(
      key,
      changeFxn,
      value as string[],
      label,
      options
    );
  }
  if (typeof data_values[0] == 'boolean') {
    return dropdownPiece(
      key,
      (v,k) => changeFxn([v]),
      value.length > 0 ? value[0] as string : null,
      label,
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

export function SingleConstraintColDefLogicPiece(
  key: string,
  changeFxn: (v: SingleSelectionDefinition, k?: string) => void, // k likely becomes meaningless so should not be used in the function
  value: SingleSelectionDefinition = emptySingleSelectionDefinition,
  label: string = 'Constraint',
  onClear: undefined | Function, 
  data_summary: DataEnvelope<any>,
  all_column_options: nestedOptionSet | string[] | undefined,
  index: number,
  totalConstraints: number,
  sorted: boolean,
  color: PropTypes.Color,
) {

  // updateFxns for inner piece calls
  // 'x' is a dummy variable needed because each piece will send a dummy key in addition the 'newDef' value.
  function updatePart(newVal, key: keyof SingleSelectionDefinition) {
    const newValue = {...value};
    newValue[key] = newVal;
    changeFxn(newValue, key);
  }
  const updateLogic = (newLogic: string | null, x?: string) => {
    updatePart(newLogic, 'logic');
  };
  const updateDef = (newDef: (string | number | null)[], x?: string) => {
    updatePart(newDef, 'def');
  };
  const updateColumn = (newCol: string, x?: string) => {
    // Reset def too!
    const newValue = {...value};
    // newValue['def'] = emptyConstraintDef;
    newValue['col'] = newCol;
    changeFxn(newValue, key);
  };
  
  // const clearDef = () => {
  //   const full = {...values};
  //   // methods (Current at index, must remain/become [emptyMethod] or lose current)
  //   let next_methods = full['methods'];
  //   next_methods.splice(index, 1);
  //   if (next_methods.length == 0) {
  //     next_methods = [emptyMethod];
  //   }
  //   // logic (Current at index-1, must remain/become [[]], lose current, or lose next if clicked for 1st def)
  //   let next_logic = full['logic'];
  //   if (next_logic.length == 1) {
  //     next_logic = [[]];
  //   } else {
  //     if (index == 0) {
  //       next_logic.splice(index, 1);
  //     }
  //     if (index > 0) {
  //       next_logic.splice(index - 1, 1);
  //     }
  //   }
  //   overallChangeFxn({methods: next_methods, logic: next_logic}, key);
  // };

  const colOptions = !!all_column_options ? all_column_options : Object.keys(data_summary);
  const chosenColValues = !!value['col'] && value['col'] in data_summary ? data_summary[value['col']] : undefined;

  const pick_column_comp = Array.isArray(colOptions) ? dropdownPiece(
    key + index + '-col',
    updateColumn,
    value['col'],
    label + ' ' + (index + 1) + ', Data Target',
    colOptions,
    sorted
  ) : nestedDropdownPiece(
    key + index,
    updateColumn,
    value['col'],
    label + ' ' + (index + 1) + ', Data Target',
    colOptions
  );
  const clear_comp = !onClear || totalConstraints <= 1 ? null : (
    <Button
      color={color}
      onClick={(x) => {
        onClear();
      }}
    >
      <DeleteIcon fontSize='small' />
    </Button>
  );
  const logic_comp = index == 0 ? null : (
    dropdownPiece(
      key + index + '-logic',
      updateLogic,
      value['logic'],
      'Combination Logic',
      ['AND', 'OR'],
      false
    )
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
        {pick_column_comp}
        <div style={{paddingLeft: 10}}>
          {!!value['col'] && SingleConstraintDefPiece(
            key + index + '-def',
            updateDef,
            value['def'],
            'Data Values',
            chosenColValues,
          )}
        </div>
      </div>
      {clear_comp}
    </div>
  );
};

export function SelectionDefinitionPiece(
  key: string,
  changeFxn: (v: SelectionDefinition, k?: string) => void, // k likely becomes meaningless so should not be used in the function
  value: SelectionDefinition = emptySelectionDefinition,
  label: string = 'Selection Criteria',
  data_summary: DataEnvelope<any>,
  all_column_options: nestedOptionSet | string[] | undefined,
  sorted: boolean,
  color: PropTypes.Color = 'primary'
) {

  function removeConstraint(index: number) {
    const newValue = [...value];
    if (index < newValue.length) {
      newValue.splice(index, 1)
      if (index==0 && newValue.length > 0) {
        // Remove the now-first constraint's logic to null
        newValue[0]['logic'] = null;
      }
    }
    changeFxn(newValue, key);
  }
  function updateConstraint(constraint: SingleSelectionDefinition, index: number) {
    const newValue = [...value];
    newValue[index] = constraint;
    changeFxn(newValue, key);
  }
  function addConstraint() {
    changeFxn([...value, emptySingleSelectionDefinition], key);
  }
  
  return (
    <div>
      <InputLabel htmlFor={`${key}-selections`} shrink>{label}</InputLabel>
      <div key={`${key}-selections`} style={{
          paddingLeft: '15px',
          paddingTop: '2px'
        }}>
        {value.map((constraint, index) => {
          return SingleConstraintColDefLogicPiece(
            key + '-constraint-' + index,
            (v, k?) => updateConstraint(v, index),
            constraint,
            'Constraint',
            () => removeConstraint(index),
            data_summary,
            all_column_options,
            index,
            value.length,
            sorted,
            color);
        })}
        <br></br>
        <Button
          color={color}
          onClick={addConstraint}
        >
          Add another condition?
        </Button>
      </div>
    </div>
  );
}

// /*
// value structure: {
//   g1: method,
//   g2: method
// }
// */
// function TwoGroupSingleConstraintDefPiece(
//   value: undefined | {g1: SelectionMethod, g2: SelectionMethod},
//   data_values: undefined | (string | number | boolean)[],
//   updateFxn: Function,
//   key: string,
//   index: number,
//   label1: string = 'Group 1',
//   label2: string = 'Group 2'
// ) {

//   function updateGroup(group, value)

//   if ( !!data_values && typeof data_values[0] == 'boolean') {
//     // Case = selecting True or False as g1, so g2 gets the other.
//     dropdownPiece(
//       key + index + index,
//       updateFxn,
//       value.length > 0 ? (value[0] as string) : undefined,
//       label,
//       ['True', 'False'],
//       false
//     );
//   }
// }



// /*
// value structure: {
//   col: string,
//   g1: method,
//   g2: method
// }
// */

// // For selection of a data column and then 2 value sets of that column (presumably mapping to two groupings that will later get).
// const DFSelectColAndTwoGroupsPiece = (
//     def: (string | number | null)[] = emptyMethod,
//     index: number,
//     data_summary: DataEnvelope<any>,
//     all_column_options: nestedOptionSet | string[],
//     updateCurrent: Function,
//     overallChangeFxn: Function,
//     key: string,
//     values: DataEnvelope<(string | number | null)[][]>,
//     sorted: boolean,
//     color: PropTypes.Color
//   ) => {
//     // updateFxns for inner piece calls
//     // 'x' is a dummy variable needed because each piece will send a dummy key in addition the 'newDef' value.
//     const updateLogic = (newLogic: string | null, x: string) => {
//       updateCurrent(null, newLogic, index);
//     };
//     const updateDefTargets = (newDef: (string | number | null)[], x: string) => {
//       // Data target inners will not have the column_name
//       const col = [def[0]];
//       updateCurrent(col.concat(...newDef), null, index);
//     };
//     const updateDefColumn = (newCol: string, x: string) => {
//       // (Resets away any data targets)
//       updateCurrent([newCol], null, index);
//     };
    
//     const clearDef = () => {
//       const full = {...values};
//       // methods (Current at index, must remain/become [emptyMethod] or lose current)
//       let next_methods = full['methods'];
//       next_methods.splice(index, 1);
//       if (next_methods.length == 0) {
//         next_methods = [emptyMethod];
//       }
//       // logic (Current at index-1, must remain/become [[]], lose current, or lose next if clicked for 1st def)
//       let next_logic = full['logic'];
//       if (next_logic.length == 1) {
//         next_logic = [[]];
//       } else {
//         if (index == 0) {
//           next_logic.splice(index, 1);
//         }
//         if (index > 0) {
//           next_logic.splice(index - 1, 1);
//         }
//       }
//       overallChangeFxn({methods: next_methods, logic: next_logic}, key);
//     };
  
//     const pick_column = Array.isArray(all_column_options) ? dropdownPiece(
//       key + index,
//       updateDefColumn,
//       def[0] as string,
//       'Condition ' + (index + 1) + ', Feature',
//       all_column_options,
//       sorted
//     ) : nestedDropdownPiece(
//       key + index,
//       updateDefColumn,
//       def[0] as string,
//       'Condition ' + (index + 1) + ', Feature',
//       all_column_options
//     );
//     const clear_comp = (
//       <Button
//         color={color}
//         onClick={(x) => {
//           clearDef();
//         }}
//       >
//         <DeleteIcon fontSize='small' />
//       </Button>
//     );
//     const logic_comp =
//       index == 0
//         ? null
//         : dropdownPiece(
//             key + index + index + index,
//             updateLogic,
//             values['logic'][index - 1][0] as string,
//             'Combination Logic',
//             ['and', 'or'],
//             false
//           );
//     const def_comp = targetSelectionComponent(
//       def,
//       data_summary,
//       updateDefTargets,
//       key,
//       index
//     );
  
//     return (
//       <div
//         key={key + index}
//         style={{
//           display: 'inline-flex'
//         }}
//       >
//         <div>
//           {logic_comp}
//           {pick_column}
//           <div style={{paddingLeft: 10}}>{def_comp}</div>
//         </div>
//         {clear_comp}
//       </div>
//     );
//   };
