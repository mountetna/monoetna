
import React, { useCallback } from 'react';
import {
  arrayLevels,
  PieceBaseInputs
} from './user_input_pieces';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import DeleteIcon from '@material-ui/icons/Delete';
import InputLabel from '@material-ui/core/InputLabel';
import { DropdownPieceRct } from './dropdown_piece';
import { NestedDropdownPieceRct } from './nested_dropdown_piece';
import { DataEnvelope, OptionSet } from '../../input_types';
import { CheckboxPieceRct } from './checkbox_piece';
import { MultiselectStringPieceRct } from './multiselect_string_piece';
import { emptyNumericConstraintDef, NumericConstraint, rangePiece, RangePieceRct } from './number_pieces';

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
export type SingleSelectionDefinition = typeof emptySingleSelectionDefinition;
export type SelectionDefinition = SingleSelectionDefinition[]
export const emptySelectionDefinition: SelectionDefinition = [emptySingleSelectionDefinition]

interface SingleConstraintDefInputs extends PieceBaseInputs<SingleSelectionConstraintDef> {
  data_values: undefined | (string | number | boolean)[]
};

function SingleConstraintDefPieceRct({
  name,
  changeFxn,
  value,
  label = 'Targets',
  data_values = undefined
}: SingleConstraintDefInputs) {

  if ( !data_values || typeof data_values[0] == 'number') {
    // Special Case = subsetting on numeric data not described in the 'data_summary'
    return <Grid item>
      <RangePieceRct
        name={name}
        changeFxn={changeFxn}
        value={value.length > 0 ? value as NumericConstraint: emptyNumericConstraintDef}
        label={label}
      />
    </Grid>;
  }
  if (typeof data_values[0] == 'string') {
    // Need to build the unique options set first.
    const options = arrayLevels(data_values) as string[];
    return <Grid item>
      <MultiselectStringPieceRct
        name={name}
        changeFxn={changeFxn as ((v: string[] | null, k?: string) => void)}
        value={value as string[]}
        label={label}
        options={options}
      />
    </Grid>;
  }
  if (typeof data_values[0] == 'boolean') {
    return <Grid item>
      <DropdownPieceRct
        name={name}
        changeFxn={(v,k) => changeFxn(v!=null? [v] : [])}
        value={value.length > 0 ? value[0] as string : null}
        label={label}
        options_in={['True', 'False']}
      />
    </Grid>;
  }
  return (
    <Grid item>
      Not yet implemented Data Taype. Please followup with someone from the Data
      Library Team!
    </Grid>
  );
}

interface SingleConstraintInputs extends PieceBaseInputs<SingleSelectionDefinition> {
  data_summary: DataEnvelope<any>,
  all_column_options?: OptionSet,
  sorted?: boolean,
  color?: 'inherit' | 'primary' | 'secondary' | 'default'
  onClear?: (() => void) | undefined,
  index: number
  totalConstraints: number
}

function SingleConstraintColDefLogicRct({
  name,
  changeFxn, // k likely becomes meaningless so should not be used in the function
  value,
  label = 'Constraint',
  onClear = undefined, 
  data_summary,
  all_column_options,
  index,
  totalConstraints,
  sorted = false,
  color = 'primary',
}: SingleConstraintInputs): React.ReactElement {

  const updatePart = useCallback((newVal, key: keyof SingleSelectionDefinition, currentValue: typeof value) => {
    const newValue = {...value};
    newValue[key] = newVal;
    changeFxn(newValue, key);
  }, [changeFxn]);
  const updateLogic = useCallback((newLogic: string | null, x: string | undefined, currentValue: typeof value) => {
    updatePart(newLogic, 'logic', currentValue);
  }, [updatePart]);
  const updateDef = useCallback((newDef: (string | number | null)[], x: string | undefined, currentValue: typeof value) => {
    updatePart(newDef, 'def', currentValue);
  }, [updatePart]);
  const updateColumn = useCallback((newCol: string | null, x: string | undefined, currentValue: typeof value) => {
    // Reset def too!
    const newValue = {...currentValue};
    newValue['def'] = emptyConstraintDef;
    newValue['col'] = newCol;
    changeFxn(newValue, name);
  }, [changeFxn]);

  const colOptions = !!all_column_options ? all_column_options : Object.keys(data_summary);
  const chosenColValues = !!value['col'] && value['col'] in data_summary ? data_summary[value['col']] : undefined;

  const pick_column_comp = <Grid item>
    {Array.isArray(colOptions) ? <DropdownPieceRct
      name={name + index + '-col'}
      changeFxn={(v,k) => {updateColumn(v, k, value)}}
      value={value['col']}
      label={label + ' ' + (index + 1) + ', Data Target'}
      options_in={colOptions}
    /> : <NestedDropdownPieceRct
      name={name + index}
      changeFxn={(v,k) => {updateColumn(v, k, value)}}
      value={value['col']}
      label={label + ' ' + (index + 1) + ', Data Target'}
      nestedOptions={colOptions}
    />}
  </Grid>
  const clear_comp = !onClear || totalConstraints <= 1 ? null : (
    <Button
      color={color}
      onClick={(x) => {
        onClear();
      }}
    >
      <DeleteIcon fontSize='small'/>
    </Button>
  );
  const logic_comp = index == 0 ? null : (
    <Grid item>
      <DropdownPieceRct
        name={name + index + '-logic'}
        changeFxn={(v,k) => {updateLogic(v, k, value)}}
        value={value['logic']}
        label={'Combination Logic'}
        options_in={['AND', 'OR']}
      />
    </Grid>
  )

  return (
    <Grid
      key={name + index}
      container
      direction='row'
    >
      <Grid item xs={10} container direction='column' alignItems='stretch'>
        {logic_comp}
        {pick_column_comp}
        <Grid item style={{paddingLeft: 10}}>
          {!!value['col'] && <SingleConstraintDefPieceRct
            name={name + index + '-def'}
            changeFxn={(v,k) => {updateDef(v, k, value)}}
            value={value['def']}
            label={'Data Values'}
            data_values={chosenColValues}
          />}
        </Grid>
      </Grid>
      {clear_comp}
    </Grid>
  );
};

interface SelectionDefinitionPieceInputs extends PieceBaseInputs<SelectionDefinition> {
  data_summary: DataEnvelope<any>,
  all_column_options?: OptionSet,
  sorted?: boolean,
  color?: 'inherit' | 'primary' | 'secondary' | 'default'
}

export function SelectionDefinitionPieceRct({
  name,
  changeFxn,
  value,
  label = '',
  data_summary,
  all_column_options,
  sorted = false,
  color = 'primary',
}: SelectionDefinitionPieceInputs): React.ReactElement | null {

  const removeConstraint = useCallback((index: number, currentValue: typeof value, name: string) => {
    const newValue = [...currentValue];
    if (index < newValue.length) {
      newValue.splice(index, 1)
      if (index==0 && newValue.length > 0) {
        // Remove the now-first constraint's logic to null
        newValue[0]['logic'] = null;
      }
    }
    changeFxn(newValue, name);
  }, [changeFxn]);
  const updateConstraint = useCallback((constraint: SingleSelectionDefinition, index: number, currentValue: typeof value, name: string) => {
    const newValue = [...currentValue];
    newValue[index] = constraint;
    changeFxn(newValue, name);
  }, [changeFxn]);
  const addConstraint = useCallback((currentValue: typeof value) => {
    changeFxn([...value, emptySingleSelectionDefinition], name);
  }, [changeFxn]);

  return (
    <div>
      <InputLabel htmlFor={`${name}-selections`} shrink>{label}</InputLabel>
      <Grid key={`${name}-selections`} container direction='column' style={{
          paddingLeft: '15px',
          paddingTop: '2px'
        }}>
        {value.map((constraint, index) => {
          return <Grid item key={`${name}-selection-${index}`}>
              <SingleConstraintColDefLogicRct
                name={name + '-constraint-' + index}
                changeFxn={(v, k?) => {updateConstraint(v, index, value, name)}}
                value={constraint}
                label={'Constraint'}
                onClear={() => {removeConstraint(index, value, name)}}
                data_summary={data_summary}
                all_column_options={all_column_options}
                index={index}
                totalConstraints={value.length}
                sorted={sorted}
                color={color}
              />
            </Grid>
          }
        )}
        <Button
          color={color}
          onClick={(v, k?) => {addConstraint(value)}}
        >
          Add another condition?
        </Button>
      </Grid>
    </div>
  );
}

export function SelectionDefinitionPiece(
  key: SelectionDefinitionPieceInputs['name'],
  changeFxn: SelectionDefinitionPieceInputs['changeFxn'], // k likely becomes meaningless so should not be used in the function
  value: SelectionDefinitionPieceInputs['value'] = emptySelectionDefinition,
  label: SelectionDefinitionPieceInputs['label'] = 'Selection Criteria',
  data_summary: SelectionDefinitionPieceInputs['data_summary'],
  all_column_options: SelectionDefinitionPieceInputs['all_column_options'],
  sorted: SelectionDefinitionPieceInputs['sorted'],
  color: SelectionDefinitionPieceInputs['color'] = 'primary'
) {
  return <SelectionDefinitionPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    data_summary={data_summary}
    all_column_options={all_column_options}
    sorted={sorted}
    color={color}
  />
}

export function OptionalSelectionDefinitionPieceRct({
  name = 'optionsl-selection',
  changeFxn,
  value = false,
  label,
  data_summary,
  all_column_options,
  sorted = false,
  color = 'primary'
}: {
  name: SelectionDefinitionPieceInputs['name'],
  changeFxn: (v: SelectionDefinition | false, k?: string) => void,
  value: SelectionDefinitionPieceInputs['value'] | false,
  label: SelectionDefinitionPieceInputs['label'],
  data_summary: SelectionDefinitionPieceInputs['data_summary'],
  all_column_options: SelectionDefinitionPieceInputs['all_column_options'],
  sorted: SelectionDefinitionPieceInputs['sorted'],
  color: SelectionDefinitionPieceInputs['color']
}): React.ReactElement {

  // Put together all subsetting options
  // let subset_options: OptionSet = Object.keys(data_summary);
  // if (additional_numeric_options!=null) {
  //   subset_options = {
  //     ...additional_numeric_options,
  //     [data_summary_options_label]: key_wrap([...subset_options])}
  // }

  const startOrClear = useCallback((doSubset: boolean, x: any) => {
    const new_full = doSubset ? emptySelectionDefinition : false;
    changeFxn(new_full, name);
  }, [changeFxn]);

  // console.log(values);
  return (
    <div key={name}>
      <CheckboxPieceRct
        name={`${name}-on-off` || 'optional-selection-on-off'}
        changeFxn={startOrClear}
        value={typeof value !== 'boolean'}
        label={label || "Do subset?"}
      />
      {typeof value !== 'boolean' && <SelectionDefinitionPieceRct
        name={`${name}-selection`}
        changeFxn={(v: SelectionDefinition, k?:string) => {changeFxn(v,name)}}
        value={value as SelectionDefinition}
        label={'Cells to keep'}
        data_summary={data_summary}
        all_column_options={all_column_options}
        sorted={sorted}
        color={color}
      />}
    </div>
  );
}

export function OptionalSelectionDefinitionPiece(
  key: SelectionDefinitionPieceInputs['name'],
  changeFxn: (v: SelectionDefinition | false, k?: string) => void,
  value: SelectionDefinitionPieceInputs['value'] | false = false,
  label: SelectionDefinitionPieceInputs['label'],
  data_summary: SelectionDefinitionPieceInputs['data_summary'],
  all_column_options: SelectionDefinitionPieceInputs['all_column_options'],
  sorted: SelectionDefinitionPieceInputs['sorted'] = false,
  color: SelectionDefinitionPieceInputs['color'] = 'primary'
) {
  return <OptionalSelectionDefinitionPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    data_summary={data_summary}
    all_column_options={all_column_options}
    sorted={sorted}
    color={color}
  />
}
