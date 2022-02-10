import React from 'react';
import {DataEnvelope} from './input_types';
import { maybeOfNullable, some, withDefault, Maybe, mapSome } from '../../../../selectors/maybe';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import MultiselectStringInput from './multiselect_string';
import { Slider } from '@material-ui/core';
import StringInput from './string';
import BooleanInput from './boolean';
import {FloatInput as EtnaFloatInput} from 'etna-js/components/inputs/numeric_input';
import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';

export function val_wrap(v: any): DataEnvelope<typeof v> {
  return {'a': v}
}

export function key_wrap(k: string[]) {
  let de: DataEnvelope<string> = {};
  for (let ind = 0; ind < k.length; ind++) {
    de[k[ind]]="0";
  }
  return de;
}

/*
"Pieces" which follow define components/elements which can be used to fill in discrete parts of a user input widget.
They are named based on types of input methods.
Overall Output Structure Requirement:
  - It is assumed that the overarching widget will produce an output consisting of a hash of key/value pairs.
Some "pieces" have additional inputs, but the first 4 are always:
  - key = the output key or name of the value element that this component should fill in. Also used as the component's key for the DOM
  - changeFxn = a function, which should be defined inside the larger ui-component in which these pieces are used.
      This funciton should take in 1) the new value and 2) the target 'key' / value elements' name, then utilize onChange to log the update.
  - value = the current value of this target element
  - label = a text label to be displayed with this component
*/

export function stringInput(
  key: string = "filler", changeFxn: Function, value: string = "filler",
  label: string = 'hello') {
    return (
      <StringInput
        key={key}
        label={label}
        data={val_wrap(value)}
        value={maybeOfNullable(value)}
        onChange={(newValue) => changeFxn(withDefault(newValue,'make'), key)}
      />
    )
  };

// export function nestableDropdownInput(
//   key: string = "filler", changeFxn: Function, value: string | null,
//   label: string, options: DataEnvelope<null>) {
//    
//     return(
//       <NestedSelectAutocompleteInput
//         key={key}
//         label={label}
//         data={{options}} 
//         value={maybeOfNullable(value)}
//         onChange={(val) => changeFxn(withDefault(val, null), key)}
//       />
//     )
//   }

export function checkboxInput(
  key: string = "filler", changeFxn: Function, value: boolean = false,
  label: string) {

    return(
      <BooleanInput
        key={key}
        label={label}
        value={maybeOfNullable(value)}
        data={val_wrap(value)}
        onChange={ value => changeFxn(withDefault(value,false), key)}
      />
    )
  }

export function dropdownInput(
  key: string = "filler", changeFxn: Function, value: string | null,
  label: string, options: string[], sorted: boolean = true) {
    
    return(
      <div key={key}>
        <Autocomplete
          disablePortal
          value={value}
          onChange={(event:any, val: string) => changeFxn(val, key)}
          id={key}
          options={options}
          renderInput={(params:any) => <TextField {...params} label={label} />}
        />
      </div>
    )
  }

export function MultiselectInput(
  key: string = "filler", changeFxn: Function, value: string[] | null,
  label: string, options: string[]) {
    
    return(
      <div key={key}>
        {label}
        <MultiselectStringInput
          data={{'0': options}}
          value={value ? some(value) : null}
          onChange={(val: Maybe<string[]>) => changeFxn(withDefault(val, null), key)}
        />
      </div>
    )
  }

export function sliderInput(
  key: string = "filler", changeFxn: Function, value: number,
  label: string, min: number = 0.1, max: number = 20) {

    return(
        <div key={key}>
          {label}
          <Slider
            value={value}
            onChange={(event, newValue) => changeFxn(newValue as number, key)}
            min={min}
            max={max}
            valueLabelDisplay="auto"
          />
        </div>
    )
  }