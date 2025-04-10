import React from 'react';
import {
  checkboxPiece,
  key_wrap,
} from './user_input_pieces';
import {PropTypes} from '@material-ui/core';
import { nestedOptionSet } from './utils';
import { emptySelectionDefinition, SelectionDefinition, SelectionDefinitionPiece } from './grouping_pieces';
import { DataEnvelope } from '../../../../ui_components';

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

export function subsetDataFramePiece(
  key: string = 'filler',
  changeFxn: (v: SelectionDefinition | false, k?: string) => void,
  value: SelectionDefinition | {} = {},
  label: string,
  data_summary_ori: DataEnvelope<any>,
  all_options: DataEnvelope<any> | undefined,
  sorted = false,
  color: PropTypes.Color = 'primary'
) {
  const data_summary = {...data_summary_ori};
  // Put together all subsetting options
  // let subset_options: string[] | nestedOptionSet = Object.keys(data_summary);
  // if (additional_numeric_options!=null) {
  //   subset_options = {
  //     ...additional_numeric_options,
  //     [data_summary_options_label]: key_wrap([...subset_options])}
  // }

  const startOrClear = (doSubset: boolean, x: any) => {
    const new_full = doSubset ? emptySelectionDefinition : false;
    console.log("Filling subetting...")
    changeFxn(new_full, key);
  };

  const base = checkboxPiece(
    key,
    startOrClear,
    Object.keys(value).length != 0,
    label
  );

  const meat = !value ? null : SelectionDefinitionPiece(
    key,
    changeFxn,
    value as SelectionDefinition,
    'Cells to keep',
    data_summary,
    all_options,
    sorted,
    color
  );

  // console.log(values);
  return (
    <div key={key}>
      {base}
      {meat}
    </div>
  );
}
