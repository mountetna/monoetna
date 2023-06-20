import React, {useEffect, useMemo} from 'react';
import {DataEnvelope} from './input_types';
import { maybeOfNullable, some, withDefault, Maybe } from '../../../../selectors/maybe';
import MultiselectStringInput from './multiselect_string';
import { InputLabel, Slider } from '@material-ui/core';
import StringInput from './string';
import BooleanInput from './boolean';
import SelectAutocompleteInput from './select_autocomplete';
import FloatInput from './float';
import { Grid, TextField } from '@material-ui/core';
import NestedSelectAutocompleteInput from './nested_select_autocomplete';
import { nestedOptionSet } from './visualizations'

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
  key: string, changeFxn: Function, value: string = '',
  label: string | undefined = undefined, minWidth: number = 150) {
    return (
      <StringInput
        key={key}
        label={label}
        value={maybeOfNullable(value)}
        data={val_wrap(value)}
        minWidth={minWidth}
        onChange={ value => changeFxn(withDefault(value,''), key)}
      />
    );
  };

export function floatPiece(
  key: string, changeFxn: Function, value: number|null = null,
  label: string | undefined = undefined, minWidth: number = 150) {
    return (
      <FloatInput
        key={key}
        label={label}
        value={some(value)}
        data={val_wrap(value)}
        minWidth={minWidth}
        onChange={ value => changeFxn(withDefault(value,null), key)}
      />
    );
  };

export function checkboxPiece(
  key: string, changeFxn: Function, value: boolean = false,
  label: string, disabled: boolean = false) {
    return(
      <BooleanInput
        key={key}
        label={label}
        value={maybeOfNullable(value)}
        data={val_wrap(value)}
        onChange={ value => changeFxn(withDefault(value,false), key)}
        disabled={disabled}
      />
    );
  }

export function dropdownPiece(
  key: string, changeFxn: Function, value: string | null = null,
  label: string | undefined, options: string[] | DataEnvelope<string>, sorted: boolean = true, minWidth: number = 200, disabled: boolean = false) {
    // options = string[] where options values are both the intended values and user-facing labels
    if (Array.isArray(options)) return(
      <SelectAutocompleteInput
        key={key}
        label={label}
        value={some(value)}
        data={val_wrap(options)}
        minWidth={minWidth}
        disabled={disabled}
        onChange={ (value) => changeFxn(withDefault(value,null), key) }
      />
    );
    // options = Object with keys = label values to show the user; values = true option values that the output / saved state should hold.
    const labels = Object.keys(options)
    const toValue = (label: string | null | undefined) => label == null? null : options[label]
    const toLabel = (value: string | null) => {
      return (value == null) ? null : (Object.keys(options).find(key => options[key] === value)) as string
    };
    return(
      <SelectAutocompleteInput
        key={key}
        label={label}
        value={some(toLabel(value))}
        data={val_wrap(labels)}
        minWidth={minWidth}
        disabled={disabled}
        onChange={ (label) => changeFxn(toValue(withDefault(label,null)), key) }
      />
    );
  }

export function nestedDropdownPiece(
  key: string, changeFxn: Function, value: string | null = null,
  label: string|undefined, options: nestedOptionSet | string[], sorted: boolean = false) {
    // sorted not implemented, but kept for compatibility with inputs that might be designed for either this or dropdownPiece and given a boolean for the sorted input there.
    if (Array.isArray(options)) {
      options = key_wrap([...options])
    }
    return <NestedSelectAutocompleteInput
      key={key}
      label={label}
      value={some(value)}
      data={val_wrap(options)}
      onChange={ (value) => changeFxn(withDefault(value,null), key) }
    />
  }

export function nestedDropdownFullPathPiece(
  key: string, changeFxn: Function, value: (string | null)[] | null = null,
  label: string|undefined, options: nestedOptionSet, sorted: boolean = true, disabled: boolean = false) {
    // options can contain multiple levels of 'categoies' which lead down to 'options'.
    // 'options' are still expected to be keys of an object
    // {categories: categories: options: null}
    if (Array.isArray(options)) {
      options = key_wrap([...options])
    }
    const value_use = value==null ? [null] : value
    return <Grid 
      key={key}
      container
      direction='column'
    >
    {value_use.map( (v, index, fullValues) => {
      const path = [...fullValues].slice(0,index)
      let options_temp = {...options}
      path.forEach((value) => {
        options_temp = value == null ? options_temp : (options_temp[value]) as DataEnvelope<DataEnvelope<null> | null>
      })
      // keys of options_temp are now the options for this level, and values are null if leaves or {options: any} for a next level
      const options_level = Object.keys(options_temp)
      return <SelectAutocompleteInput
        key={key+index}
        label={index==0 ? label : undefined}
        value={some(v)}
        data={val_wrap(options_level)}
        disabled={disabled}
        onChange={ (value) => {
          // Trim to the level where the new selection was made
          const val = withDefault(value,null)
          path.push(val)
          // Add new level if selection was not a leaf
          if (val != null && options_temp[val]!=null) {path.push(null)}
          changeFxn(path, key) }
        }
      />
    })}
    </Grid>
  }

export function MultiselectAfterDataChoicePiece(
  key: string, changeFxn: Function, value: string[] | null = null,
  label: string,
  full_data: DataEnvelope<any[]>,
  data_target: string | null, // The name of a column/key of full_data which the user has chosen as the target of this ui piece, or null if not chosen yet.
  data_target_label: string, // The label of the ui-piece where the user selects data_target 
  discrete_data: string[] | nestedOptionSet
) {

  const canReorder = data_target != null && discrete_data.includes(data_target) && full_data != null
  const levels = canReorder ? arrayLevels(Object.values(full_data[data_target as string])) : null
  return(
    levels != null ? multiselectPiece(key, changeFxn, value, label, levels) :
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

export function multiselectPiece(
  key: string, changeFxn: Function, value: string[] | null = null,
  label: string, options: string[]) {
    useEffect(() => {
      if (options != null && value != null && value.length != 0) {
        // Needed if given options change and chosen ones aren't captured
        let needs_reset = 
            value.filter(
              (val) => !options.includes(val)
            ).length > 0;
        if (needs_reset) changeFxn(withDefault([], null), key);
      }
    }, [options]);
    
    return(
      <div key={key} style={{paddingTop: 8}}>
        <InputLabel htmlFor={'multiselect-'+key} shrink>{label}</InputLabel>
        <MultiselectStringInput
          key={'multiselect-'+key}
          data={{'0': options}}
          value={maybeOfNullable(value)}
          onClear={() => changeFxn([], key)}
          onChange={(val: Maybe<string[]>) => changeFxn(withDefault(val, null), key)}
        />
      </div>
    );
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

export function rangePiece(
  key: string, changeFxn: Function, value: (string|number|null)[] = ['exactly', null, 'below', null],
  label: string) {
    const updateSlot = (newValue: string|number|null, slot: number, current_full = value) => {
      let next_full = [...current_full];
      next_full[slot] = newValue;
      return next_full;
    };
    
    return(
      <div key={key}>
        <div style={{display: 'inline-flex'}}>
          {dropdownPiece(
            key+'_lower_bound_type', (newValue: string | null) => changeFxn(updateSlot(newValue, 0), key), value[0] as string,
            label + ', From', ['exactly','above'], true, 120)}
          {floatPiece(
            key+'_lower_value', (newValue: number | null) => changeFxn(updateSlot(newValue, 1), key), value[1] as number,
            'Min-value', 120)}
        </div>
        <div style={{display: 'inline-flex'}}>
          {dropdownPiece(
            key+'_upper_bound_type', (newValue: string | null) => changeFxn(updateSlot(newValue, 2), key), value[2] as string,
            'To', ['exactly','below'], true, 120)}
          {floatPiece(
            key+'_upper_value', (newValue: number | null) => changeFxn(updateSlot(newValue, 3), key), value[3] as number,
            'Max-value', 120)}
        </div>
      </div>
    );
  }

export function reductionSetupPiece(
  key: string, changeFxn: Function, value: (string|null)[] = [null, '1', '2'],
  label: string[] = ['Dimensionality Reduction (DR)', 'x-axis DR Compenent', 'y-axis DR Component'],
  reduction_opts: DataEnvelope<string[]>) {
  
  if (reduction_opts == null) return null

  const disable_dims = value[0]==null
  function changeReduction(newElement: string|null) {
    return [newElement, '1', '2']
  }
  function changeDim(newElement: string|null, dim: 1|2) {
    let newValue = [...value]
    newValue[dim] = newElement
    return newValue
  }
  // console.log({reduction_opts})
  return(
    <Grid 
      key={key}
      container
      direction='column'
    >
      <Grid item>
        {dropdownPiece(
          key+'-reduction', (newElement: string | null) => changeFxn(changeReduction(newElement), key), value[0],
          label[0], Object.keys(reduction_opts), false, 200, false)}
      </Grid>
      <Grid item>
        {dropdownPiece(
          key+'-dimx', (newElement: string | null) => changeFxn(changeDim(newElement, 1), key), value[1],
          label[1], value[0]==null ? ['1', '2'] : reduction_opts[value[0]], false, 200, disable_dims)}
      </Grid>
      <Grid item>
        {dropdownPiece(
          key+'-dimy', (newElement: string | null) => changeFxn(changeDim(newElement, 2), key), value[2],
          label[2], value[0]==null ? ['1', '2'] : reduction_opts[value[0]], false, 200, disable_dims)}
      </Grid>
    </Grid>
  )

}
