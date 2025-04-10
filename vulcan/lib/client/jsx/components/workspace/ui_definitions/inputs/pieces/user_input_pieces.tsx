import React, { useEffect } from 'react';
import {DataEnvelope} from '../../input_types';
import { maybeOfNullable, some, withDefault, Maybe } from '../../../../../selectors/maybe';
import InputLabel from '@material-ui/core/InputLabel';
import Slider from '@material-ui/core/Slider';
import FloatInput from '../components/float';
import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';
import NestedDropdownInput from '../components/nested_dropdown';
import { nestedOptionSet } from './utils'
import NestedDropdownMultiChoicePieceRct from './nested_dropdown_multi_choice_piece';
import { FloatPiece } from './number_pieces';

export function val_wrap(v: any): DataEnvelope<typeof v> {
  return {'a': v};
}

export function key_wrap(k: string[]) {
  let de: DataEnvelope<null> = {};
  for (let ind = 0; ind < k.length; ind++) {
    de[k[ind]]=null;
  }
  return de;
}

export function arrayLevels(original: any[]) {
  function onlyUnique(value: any, index: number, self: any) {
    return self.indexOf(value) === index;
  }
  return Array.from(original).filter(onlyUnique);
}

/*
"Pieces" which follow define components/elements which can be used to fill in discrete parts of a user input widget.
They are named based on types of input methods.
Overall Output Structure Requirement:
  - It is assumed that the overarching widget will produce an output consisting of a hash of key/value pairs.
Some "pieces" have additional inputs, but the first 4 are ALWAYS:
  - key = the output key or name of the value element that this component should fill in. Also used as the component's key for the DOM
  - changeFxn = a function, which should be defined inside the larger ui-component in which these pieces are used.
      This funciton should take in 1) the new value and 2) the target 'key' / value elements' name, then utilize onChange to log the update.
  - value = the current value of this target element
  - label = a text label to be displayed with this component
*/

export function stringPiece(
  key: string, changeFxn: (v: any, k?:string) => void, value: string = '',
  label: string | undefined = undefined, minWidth: number = 150) {
    return (
      <div key={key} style={{paddingTop: label ? 8 : 0}}>
        <TextField
          value={value}
          multiline
          label={label}
          InputLabelProps={{ shrink: true }}
          onChange={(event) => changeFxn(event.target.value, key)}
          size="small"
          style={{minWidth: minWidth || 200}}
        />
      </div>
    );
  };

export function checkboxPiece(
  key: string, changeFxn: (v: boolean, k?:string) => void, value: boolean = false,
  label: string,
  disabled: boolean = false,
  labelPlacement: FormControlLabelProps['labelPlacement'] ='end') {

  return <CheckboxPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    disabled={disabled}
    labelPlacement={labelPlacement}
  />
};

export function CheckboxPieceRct({
  name, changeFxn, value = false,
  label,
  disabled = false,
  labelPlacement = 'end'}: {
    name: string,
    changeFxn: (v: boolean, k?:string) => void,
    value: boolean,
    label: string,
    disabled?: boolean,
    labelPlacement?: FormControlLabelProps['labelPlacement']
  }): React.ReactElement {

  const inner = <Checkbox
    onChange={(e: any) => changeFxn(!value, name)}
    checked={value}
    inputProps={{ 'aria-label': 'controlled' }}
    disabled={disabled}
  />;
  if (label) {
    return (
      <FormControlLabel key={name} control={inner} label={label} labelPlacement={labelPlacement}/>
    );
  }

  return inner;
};

import ListInput from 'etna-js/components/inputs/list_input';
import DropdownAutocompleteInput from 'etna-js/components/inputs/dropdown_autocomplete_wrapper';
import DropdownPiece from './dropdown_piece';
import Checkbox from '@material-ui/core/Checkbox';
import FormControlLabel, { FormControlLabelProps } from '@material-ui/core/FormControlLabel';
export function MultiselectPiece(
  key: string, changeFxn: Function, value: string[] | null = null,
  label: string, options: string[],
  onAll: string[] = options,
  onClear: string[] = [],
  maxItems: number = 25
) {

  if (value != null && value.length > 0) {
    if (value.filter((val) => !(val==='' || options.includes(val))).length > 0) changeFxn(null, key);
  }

  return(
    <div key={key} style={{paddingTop: 8}}>
      <InputLabel htmlFor={'multiselect-'+key} shrink>{label}</InputLabel>
      <ListInput
        key={'multiselect-'+key}
        placeholder='Select items from the list'
        className='link_text'
        values={value!=null ? value : []}
        itemInput={DropdownAutocompleteInput}
        list={options}
        onChange={(val: string[]) => changeFxn(val, key)}
        onAll={() => changeFxn(onAll, key)}
        onClear={() => changeFxn(onClear, key)}
        maxItems={maxItems}
      />
    </div>
  );
}

  export function MultiselectAfterDataChoicePiece(
    key: string, changeFxn: Function, value: string[] | null = null,
    label: string,
    full_data: DataEnvelope<any[]>,
    data_target: string | null, // The name of a column/key of full_data which the user has chosen as the target of this ui piece, or null if not chosen yet.
    data_target_label: string, // The label of the ui-piece where the user selects data_target 
    discrete_data: string[]
  ) {
  
    const canReorder = data_target != null && discrete_data.includes(data_target) && full_data != null
    const levels = canReorder ? arrayLevels(Object.values(full_data[data_target as string])) : null
    return(
      levels != null ? MultiselectPiece(key, changeFxn, value, label, levels) :
        <div key={key} style={{paddingTop: 8}}>
          <TextField
            key={'multiselect-'+key}
            label={label}
            InputLabelProps={{ shrink: true }}
            size="small"
            disabled
            fullWidth
            placeholder={'Awaiting '+data_target_label+' choice'}
          />
        </div>
    )
  }

export function sliderPiece(
  key: string, changeFxn: Function, value: number,
  label: string, min: number = 0.1, max: number = 20, stepSize: number | undefined = undefined) {

    return(
        <div key={key} style={{paddingTop: 8}}>
          <InputLabel htmlFor={'slider-'+key} shrink>{label}</InputLabel>
          <Slider
            key={'slider-'+key}
            value={value}
            onChange={(event, newValue) => changeFxn(newValue as number, key)}
            min={min}
            max={max}
            step={stepSize}
            valueLabelDisplay="auto"
          />
        </div>
    );
  }

export type NumericConstraint = ['exactly' | 'above', number | null, 'exactly' | 'below', number | null];
export const emptyNumericConstraintDef: NumericConstraint = ['exactly', 0, 'below', 0];
export function rangePiece(
  key: string, changeFxn: Function, value: NumericConstraint,
  label: string) {
    const updateSlot = (newValue: 'exactly' | 'above' | 'below' | number | null, slot: number, current_full = value) => {
      let next_full = [...current_full];
      next_full[slot] = newValue;
      return next_full;
    };
    
    return(
      <div>
        {label==='' ? null : <InputLabel htmlFor={`${key}-range`} shrink>{label}</InputLabel>}
        <div key={`${key}-range`}>
          <div style={{display: 'inline-flex'}}>
            {DropdownPiece(
              key+'_lower_bound_type', (newValue: string | null) => changeFxn(updateSlot(newValue as 'exactly' | 'above', 0), key), value[0] as string,
              'From', ['exactly','above'], true, 120)}
            {FloatPiece(
              key+'_lower_value', (newValue: number | null) => changeFxn(updateSlot(newValue as number | null, 1), key), value[1] as number,
              'Min-value', 120)}
          </div>
          <div style={{display: 'inline-flex'}}>
            {DropdownPiece(
              key+'_upper_bound_type', (newValue: string | null) => changeFxn(updateSlot(newValue as 'exactly' | 'below', 2), key), value[2] as string,
              'To', ['exactly','below'], true, 120)}
            {FloatPiece(
              key+'_upper_value', (newValue: number | null) => changeFxn(updateSlot(newValue as number | null, 3), key), value[3] as number,
              'Max-value', 120)}
          </div>
        </div>
      </div>
    );
  }

export function ReductionSetupPiece(
  key: string, changeFxn: Function, value: (string|null)[] = [null, '1', '2'],
  label: string[] = ['Dimensionality Reduction (DR)', 'x-axis DR Compenent', 'y-axis DR Component'],
  reduction_opts: DataEnvelope<string[]>) {
  
  function changeReduction(newElement: string|null) {
    return [newElement, '1', '2']
  }
  function changeDim(newElement: string|null, dim: 1|2) {
    let newValue = [...value]
    newValue[dim] = newElement
    return newValue
  }

  useEffect(()=>{
    // Default to _Recommended_ (which comes from 'umap_name' attribute of the dataset record)
    if (value[0]==null && reduction_opts != null && Object.keys(reduction_opts).includes('_Recommended_')) {
      changeReduction('_Recommended_')
    }
  }, [value, reduction_opts])

  const error_with_input = !reduction_opts;

  const disable_dims = value[0]==null || error_with_input;
  // console.log({reduction_opts})
  return(
    <Grid 
      key={key}
      container
      direction='column'
    >
      <Grid item>
        {error_with_input ?
          DropdownPiece(
            key+'-reduction', (newElement: string | null) => changeFxn(changeReduction(newElement), key), 'Error in input setup',
            label[0], ['Error in input setup'], false, 200, true) :
          DropdownPiece(
            key+'-reduction', (newElement: string | null) => changeFxn(changeReduction(newElement), key), value[0],
            label[0], Object.keys(reduction_opts), false, 200, false)}
      </Grid>
      <Grid item style={{paddingLeft: 10}}>
        {DropdownPiece(
          key+'-dimx', (newElement: string | null) => changeFxn(changeDim(newElement, 1), key), value[1],
          label[1], value[0]==null ? ['1', '2'] : reduction_opts[value[0]] || ['1'], false, 200, disable_dims)}
      </Grid>
      <Grid item style={{paddingLeft: 10}}>
        {DropdownPiece(
          key+'-dimy', (newElement: string | null) => changeFxn(changeDim(newElement, 2), key), value[2],
          label[2], value[0]==null ? ['1', '2'] : reduction_opts[value[0]] || ['2'], false, 200, disable_dims)}
      </Grid>
    </Grid>
  )

}
