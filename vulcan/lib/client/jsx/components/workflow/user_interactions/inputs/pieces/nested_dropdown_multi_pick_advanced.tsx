import React, {useEffect, useMemo, useState, useCallback} from 'react';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import AddBoxIcon from '@material-ui/icons/AddBox';
import Tooltip from '@material-ui/core/Tooltip';
import { nestedOptionSet } from '../visualizations'
import { ReorderCollapsiblePiece } from '../reorder_piece';
import { nestedDropdownMultiPickPiece } from '../user_input_pieces';
import SelectionsFromPasteModal, { fuseSearchSets } from './select_from_paste_modal';
import { flattenOptionPaths, leafParentPaths, pathValues, sep } from '../nested_select_autocomplete_multi_choice';

function val_wrap_in_val_key(v: string[]): {val: string}[] {
  return v.map((value) => {
    return {val: value};
  });
}

export default function NestedDropdownMultiPickAdvanced(
  key: string, changeFxn: Function, value: string[] = [] as string[],
  label: string|undefined, options: nestedOptionSet | string[], sorted: boolean = false) {
  // sorted not implemented, but kept for compatibility with inputs that might be designed for either this or dropdownPiece and given a boolean for the sorted input there.

  const [open, setOpen] = useState(false);
  const [loading, setLoading] = useState(false);

  // Output function for 'Bulk Add' modal, merge with current values
  const mergeValues = useCallback( (newValues: string[]) => {
    if (value==null) {
      changeFxn(newValues, key);
    } else {
      changeFxn(value.concat(newValues.filter((val)=>!value.includes(val))), key);
    }
  }, [value, key, changeFxn]);

  // Convert the given options to [{val: option1}, {val: option2}, ...] per each 'option set'
  const optionSets: fuseSearchSets = useMemo(()=>{
    if (Array.isArray(options)) return {Options: val_wrap_in_val_key(options)};
    const sets = leafParentPaths(flattenOptionPaths(options), sep);
    const output: fuseSearchSets = {};
    for (const set of sets) {
      output[set] = val_wrap_in_val_key(pathValues(set, options, sep));
    }
    return output
  }, [options])
  
  return <Grid container direction='column'>
    <Grid item container alignItems='center'>
      <Grid item xs={10}>
        {nestedDropdownMultiPickPiece(`${key}-selection`, (value: string[]) => changeFxn(value, key), value, `${label} - selection`, options)}
      </Grid>
      <Grid item>
        <Tooltip title="Bulk Add">
          <IconButton
            aria-label="Bulk Add"
            onClick={() => {
              setOpen(true);
              setLoading(true);
            }}
            color='secondary'
            disabled={loading}
            variant='text'
          >
            <AddBoxIcon/>
          </IconButton>
        </Tooltip>
      </Grid>
    </Grid>
    <Grid item style={{paddingLeft: '12px'}}>
      {ReorderCollapsiblePiece(`${key}-reorder`, (value: string[]) => changeFxn(value, key), value, 'Reorder selections?')}
    </Grid>
    <Dialog
      open={open}
      onClose={(event: object, reason: string) => {if (reason!="backdropClick") setOpen(false)}}
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
        onClose={() => {setOpen(false)}}
        genomicData={true}
      />
    </Dialog>
  </Grid>
}