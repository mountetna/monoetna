import React from 'react';
import {DataEnvelope} from './input_types';
import { maybeOfNullable, some, withDefault, Maybe } from '../../../../selectors/maybe';
import MultiselectStringInput from './multiselect_string';
import { Slider } from '@material-ui/core';
import StringInput from './string';
import BooleanInput from './boolean';
import SelectAutocompleteInput from './select_autocomplete';
import FloatInput from './float';

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
Some "pieces" have additional inputs, but the first 4 are ALWAYS:
  - key = the output key or name of the value element that this component should fill in. Also used as the component's key for the DOM
  - changeFxn = a function, which should be defined inside the larger ui-component in which these pieces are used.
      This funciton should take in 1) the new value and 2) the target 'key' / value elements' name, then utilize onChange to log the update.
  - value = the current value of this target element
  - label = a text label to be displayed with this component
*/

export function stringPiece(
  key: string, changeFxn: Function, value: string = "",
  label: string | undefined = undefined, minWidth: number = 150) {
    return (
      <StringInput
        key={key}
        label={label}
        value={maybeOfNullable(value)}
        data={val_wrap(value)}
        minWidth={minWidth}
        onChange={ value => changeFxn(withDefault(value,""), key)}
      />
    )
  };

export function floatPiece(
  key: string, changeFxn: Function, value: number|null = null,
  label: string | undefined = undefined, minWidth: number = 150) {
    return (
      <FloatInput
        key={key}
        label={label}
        value={maybeOfNullable(value)}
        data={val_wrap(value)}
        minWidth={minWidth}
        onChange={ value => changeFxn(withDefault(value,null), key)}
      />
    )
  };

export function checkboxPiece(
  key: string, changeFxn: Function, value: boolean = false,
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

export function dropdownPiece(
  key: string, changeFxn: Function, value: string | null = null,
  label: string|undefined, options: string[], sorted: boolean = true, minWidth: number = 200) {
    return(
      <SelectAutocompleteInput
        key={key}
        label={label}
        value={some(value)}
        data={val_wrap(options)}
        minWidth={minWidth}
        onChange={ (value) => changeFxn(withDefault(value,null), key) }
      />
    )
  }

export function MultiselectPiece(
  key: string, changeFxn: Function, value: string[] | null = null,
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

export function sliderPiece(
  key: string, changeFxn: Function, value: number,
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

export function rangePiece(
  key: string, changeFxn: Function, value: (string|number|null)[] = ["exactly", null, "below", null],
  label: string) {
    const updateSlot = (newValue: string|number|null, slot: number, current_full = value) => {
      let next_full = [...current_full]
      next_full[slot] = newValue
      return next_full
    }
    
    return(
      <div key={key}>
        {label}
        <div style={{display: 'inline-flex', paddingTop:8}}>
          From:
          {dropdownPiece(
            key+'_lower_bound_type', (newValue: string | null) => changeFxn(updateSlot(newValue, 0), key), value[0] as string,
            undefined, ["exactly","above"], true, 120)}
          {floatPiece(
            key+'_lower_value', (newValue: number | null) => changeFxn(updateSlot(newValue, 1), key), value[1] as number,
            'Min-value', 120)}
        </div>
        <div style={{display: 'inline-flex', paddingTop:8}}>
          To:
          {dropdownPiece(
            key+'_upper_bound_type', (newValue: string | null) => changeFxn(updateSlot(newValue, 2), key), value[2] as string,
            undefined, ["exactly","below"], true, 120)}
          {floatPiece(
            key+'_upper_value', (newValue: number | null) => changeFxn(updateSlot(newValue, 3), key), value[3] as number,
            'Max-value', 120)}
        </div>
      </div>
    )
  }
