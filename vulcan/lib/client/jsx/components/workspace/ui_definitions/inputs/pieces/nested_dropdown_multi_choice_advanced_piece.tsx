import React, {useEffect, useMemo, useState, useCallback} from 'react';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import LibraryAddRoundedIcon from '@material-ui/icons/LibraryAddRounded';
import Tooltip from '@material-ui/core/Tooltip';
import { ReorderCollapsiblePieceRct } from './reorder_piece';
import SelectionsFromPasteModal, { fuseSearchSets } from './select_from_paste_modal';
import NestedDropdownMultiChoicePieceRct, { flattenOptionPaths, leafParentPaths, NestedDropdownMultiChoicePieceInputs, pathValues, sep } from './nested_dropdown_multi_choice_piece';
import { key_wrap } from './user_input_pieces';
import { OptionSet } from '../../input_types';

function val_wrap_in_val_key(v: string[]): {val: string}[] {
  return v.map((value) => {
    return {val: value};
  });
}

const buttonStyles = makeStyles((theme) => ({
  gridFullHieght: {
    display: 'flex',
    height: '100%',
  },
  centerHeightChildren: {
    display: 'flex',
    alignItems: 'center'
  }
}));

interface NestedDropdownMultiChoiceAdvancedInputs extends NestedDropdownMultiChoicePieceInputs {
  letReorder: boolean,
  letBulkAdd: boolean,
};

export default function NestedDropdownMultiChoiceAdvancedPiece(
  key: NestedDropdownMultiChoiceAdvancedInputs['name'],
  changeFxn: NestedDropdownMultiChoiceAdvancedInputs['changeFxn'],
  value: NestedDropdownMultiChoiceAdvancedInputs['value'] = [] as string[],
  label: NestedDropdownMultiChoiceAdvancedInputs['label'],
  options: OptionSet,
  letReorder: NestedDropdownMultiChoiceAdvancedInputs['letReorder'] = true,
  letBulkAdd: NestedDropdownMultiChoiceAdvancedInputs['letBulkAdd'] = true,
  testId: NestedDropdownMultiChoiceAdvancedInputs['testId'],
  sorted: NestedDropdownMultiChoiceAdvancedInputs['sorted']
): React.ReactElement {
  return <NestedDropdownMultiChoiceAdvancedPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    options_in={options}
    letReorder={letReorder}
    letBulkAdd={letBulkAdd}
    testId={testId}
    sorted={sorted}
  />
}

export function NestedDropdownMultiChoiceAdvancedPieceRct({
  name,
  changeFxn,
  value = [] as string[],
  label,
  options_in,
  letReorder = true,
  letBulkAdd = true,
  testId = undefined,
  sorted = true
}: NestedDropdownMultiChoiceAdvancedInputs) {

  const classes = buttonStyles();

  const [open, setOpen] = useState(false);
  const [loading, setLoading] = useState(false);

  // Output function for 'Bulk Add' modal, merge with current values
  const mergeValues = useCallback( (newValues: string[]) => {
    if (value==null) {
      changeFxn(newValues, name);
    } else {
      changeFxn(value.concat(newValues.filter((val)=>!value.includes(val))), name);
    };
  }, [value, name, changeFxn]);

  // Convert the given options to format needed for fuse.search(), [{val: option1}, {val: option2}, ...] per each 'option set'
  const optionSets: fuseSearchSets = useMemo(()=>{
    if (Array.isArray(options_in)) return {Options: val_wrap_in_val_key(options_in)};
    const sets = leafParentPaths(flattenOptionPaths(options_in), sep);
    const output: fuseSearchSets = {};
    for (const set of sets) {
      output[set] = val_wrap_in_val_key(pathValues(set, options_in, sep));
    };
    return output;
  }, [options_in]);

  let bulkAdd: any = null;
  if (letBulkAdd) {

    bulkAdd = <Grid item className={classes.gridFullHieght}>
      <Tooltip title="Bulk Add" placement='right'>
        <IconButton
          aria-label="Bulk Add"
          onClick={() => {
            setOpen(true);
            setLoading(true);
          }}
          color='secondary'
          disabled={loading}
          variant='outlined'
          className={classes.centerHeightChildren}
        >
          <LibraryAddRoundedIcon style={{transform: 'rotate(90deg)'}}/>
        </IconButton>
      </Tooltip>
      <Dialog
        open={open}
        onClose={(event: object, reason: string) => {if (reason!='backdropClick') setOpen(false);}}
        maxWidth='xl'
        TransitionProps={{
          onEntered: () => {
            setLoading(false);
          }
        }}
        disableAutoFocus={true}
        disableEnforceFocus={true}
      >
        <SelectionsFromPasteModal
          options={optionSets}
          onSave={mergeValues}
          onClose={() => {setOpen(false);}}
          genomicData={true}
        />
      </Dialog>
    </Grid>;
  };

  return <Grid container direction='column'>
    <Grid item container alignItems='center'>
      <Grid item xs={letBulkAdd ? 10 : false}>
        <NestedDropdownMultiChoicePieceRct
          name={`${name}-selection`}
          changeFxn={(v, k?) => {changeFxn(v, name)}}
          value={value}
          label={`${label} - selection`}
          options_in={options_in}
          testId={testId}
        />
      </Grid>
      {bulkAdd}
    </Grid>
    {letReorder && <Grid item style={{paddingLeft: '12px'}}>
      <ReorderCollapsiblePieceRct
        name={`${name}-reorder`}
        changeFxn={(value: string[], k?) => changeFxn(value, name)}
        value={!!value ? value : []}
        label='Reorder selections?'
      />
    </Grid>}
  </Grid>;
};
