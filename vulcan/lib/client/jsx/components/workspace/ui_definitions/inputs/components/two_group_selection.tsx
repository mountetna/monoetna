import React from 'react';
import {DataEnvelope, OptionSet, WithInputParams} from '../../input_types';
import {some} from '../../../../../selectors/maybe';
import {useSetsDefault} from '../../useSetsDefault';
import { emptySelectionDefinition, SelectionDefinitionPiece, SelectionDefinition, SelectionDefinitionPieceRct } from '../pieces/grouping_pieces';
import Grid from '@material-ui/core/Grid';
import InputLabel from '@material-ui/core/InputLabel';

export default function TwoGroupSelection({
  data,
  label,
  onChange,
  ...props
}: WithInputParams<{}, DataEnvelope<any>, {data_summary: DataEnvelope<any>, all_column_options?: OptionSet}>) {
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

  const key = label || 'two-group-selection';

  // console.log({g1_value})
  // console.log({g2_value})

  return <div>
    {label!=null && <InputLabel htmlFor={key} shrink>{label}</InputLabel>}
    <Grid key={key} container direction='column'>
      <Grid item>
        {!!g1_value && <SelectionDefinitionPieceRct
          name='g1'
          changeFxn={(v: SelectionDefinition) => onChange({ g1: some(v), g2: some(g2_value) })}
          value={g1_value}
          label='Group 1 Selection Criteria'
          data_summary={data['data_summary']}
          all_column_options={'all_column_options' in data ? data['all_column_options'] : undefined}
          color='primary'
        />}
      </Grid>
      <Grid item>
        {!!g2_value && <SelectionDefinitionPieceRct
          name='g2'
          changeFxn={(v: SelectionDefinition) => onChange({ g1: some(g1_value), g2: some(v) })}
          value={g2_value}
          label='Group 2 Selection Criteria'
          data_summary={data['data_summary']}
          all_column_options={'all_column_options' in data ? data['all_column_options'] : undefined}
          color='primary'
        />}
      </Grid>
    </Grid>
  </div>
}