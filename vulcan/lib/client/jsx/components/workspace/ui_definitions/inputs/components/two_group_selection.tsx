import React from 'react';
import {DataEnvelope, WithInputParams} from '../../input_types';
import {some} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { emptySelectionDefinition, SelectionDefinitionPiece, SelectionDefinition } from '../pieces/grouping_pieces';
import { nestedOptionSet } from '../pieces/utils';
import Grid from '@material-ui/core/Grid';

export default function TwoGroupSelection({
  data,
  onChange,
  ...props
}: WithInputParams<{}, DataEnvelope<any>, {data_summary: DataEnvelope<any>, all_column_options?: nestedOptionSet | string[]}>) {
  const g1_value = useSetsDefault(
    emptySelectionDefinition,
    props.value,
    onChange,
    'g1'
  );
  const g2_value = useSetsDefault(
    emptySelectionDefinition,
    props.value,
    onChange,
    'g2'
  );

  // console.log({g1_value})
  // console.log({g2_value})

  return (
    <Grid container direction='column'>
      <Grid item>
        {!!g1_value && SelectionDefinitionPiece(
          'g1',
          (v: SelectionDefinition) => onChange({g1: some(v), g2: some(g2_value)}),
          g1_value,
          'Group 1 Selection Criteria',
          data['data_summary'],
          'all_column_options' in data ? data['all_column_options'] : undefined,
          false,
          'primary'
        )}
      </Grid>
      <Grid item>
        {!!g2_value && SelectionDefinitionPiece(
          'g2',
          (v: SelectionDefinition) => onChange({g1: some(g1_value), g2: some(v)}),
          g2_value,
          'Group 2 Selection Criteria',
          data['data_summary'],
          'all_column_options' in data ? data['all_column_options'] : undefined,
          false,
          'primary'
        )}
      </Grid>
    </Grid>
  )
}