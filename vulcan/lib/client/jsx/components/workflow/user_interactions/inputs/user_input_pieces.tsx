import React from 'react';
import {DataEnvelope} from './input_types';
import { maybeOfNullable, some, withDefault, Maybe } from '../../../../selectors/maybe';
import MultiselectStringInput from './multiselect_string';
import { InputLabel, Paper, Slider } from '@material-ui/core';
import StringInput from './string';
import BooleanInput from './boolean';
import SelectAutocompleteInput from './select_autocomplete';
import FloatInput from './float';
import {DragDropContext, Droppable, Draggable} from 'react-beautiful-dnd';
import { arrayLevels } from './subsetDataFrame_piece';

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
        value={some(value)}
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

export function multiselectPiece(
  key: string, changeFxn: Function, value: string[] | null = null,
  label: string, options: string[]) {
    
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
    )
  }

export function sliderPiece(
  key: string, changeFxn: Function, value: number,
  label: string, min: number = 0.1, max: number = 20) {

    return(
        <div key={key} style={{paddingTop: 8}}>
          <InputLabel htmlFor={'slider-'+key} shrink>{label}</InputLabel>
          <Slider
            key={'slider-'+key}
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
        <div style={{display: 'inline-flex'}}>
          {dropdownPiece(
            key+'_lower_bound_type', (newValue: string | null) => changeFxn(updateSlot(newValue, 0), key), value[0] as string,
            label + ", From", ["exactly","above"], true, 120)}
          {floatPiece(
            key+'_lower_value', (newValue: number | null) => changeFxn(updateSlot(newValue, 1), key), value[1] as number,
            'Min-value', 120)}
        </div>
        <div style={{display: 'inline-flex'}}>
          {dropdownPiece(
            key+'_upper_bound_type', (newValue: string | null) => changeFxn(updateSlot(newValue, 2), key), value[2] as string,
            "To", ["exactly","below"], true, 120)}
          {floatPiece(
            key+'_upper_value', (newValue: number | null) => changeFxn(updateSlot(newValue, 3), key), value[3] as number,
            'Max-value', 120)}
        </div>
      </div>
    )
  }

export function reorderPiece(
  key: string, changeFxn: Function, value: string | string[] = "unordered",
  label: string, full_data: DataEnvelope<any[]>, data_target: string | null) {

  let chosen_case = Array.isArray(value) ? "custom" : value

  const case_dropdown = <SelectAutocompleteInput
      key={key}
      label={label}
      disabled={data_target == null}
      value={some(chosen_case)}
      data={{a: ['unordered', 'increasing', 'decreasing', 'custom']}}
      onChange={() => {}}
      onChangeOverride={ (e, picked) => {
        console.log({picked})
        const newValue = picked=='custom' ? arrayLevels(Object.values(full_data[data_target as string])) : picked
        changeFxn(newValue, key)
      }}
    />
  
  const LevelComponent = (props: any) => {
    return(
      <Draggable
        draggableId={`level-${props.levelIndex}`}
        index={props.levelIndex}
      >
        {(provided: any) => (
          <div
            ref={provided.innerRef}
            {...provided.draggableProps}
            {...provided.dragHandleProps}
          >
            {props.level}
          </div>
        )}
      </Draggable>
    )
  }

  const handleOnDragEnd = (result: any) => {
    if (!result.destination) {
      return;
    }

    if (result.destination.index === result.source.index) {
      return;
    }

    const newValues = Array.from(value)
    const [removed] = newValues.splice(result.source.index, 1);
    newValues.splice(result.destination.index, 0, removed);

    changeFxn(newValues, key)
  }
  
  const reorder_custom = !Array.isArray(value) ? null :
    <DragDropContext onDragEnd={handleOnDragEnd}>
      <Droppable droppableId='columns'>
        {(provided: any) => (
          <Paper ref={provided.innerRef} {...provided.droppableProps}>
            {(value as string[]).map((level: string, index: number) => {
              return (
                <LevelComponent
                  key={index}
                  level={level}
                  levelIndex={index}
                />
              );
            })}
            {provided.placeholder}
          </Paper>
        )}
      </Droppable>
    </DragDropContext>
  
  console.log({value})
  
  return(
    <div key={key} style={{paddingTop: 8}}>
      {case_dropdown}
      {reorder_custom}
    </div>
  )
}