import React, { useCallback, useEffect, useState } from 'react';
import TextField from '@material-ui/core/TextField';
import { PieceBaseInputs } from './user_input_pieces';
import InputLabel from '@material-ui/core/InputLabel';
import Slider from '@material-ui/core/Slider';
import { DropdownPieceRct } from './dropdown_piece';
import Grid from '@material-ui/core/Grid';

function parseIntBetter(s: string) {
  const parsed = parseInt(s, 10);
  return parsed + '' == s ? parsed : NaN;
}

function from_num(num: number | null, float = false) {
  if (float) {
    return num == null ?
      '0.0' :
      num.toString().includes('.') ?
        num.toString() :
        num.toString() + '.0';
  }
  return num == null ? '' : num.toString();
}

interface NumberPieceInputs extends PieceBaseInputs<number | null> {
  minWidth?: number;
  integer: boolean;
}

export function NumberPieceRct({
  name,
  changeFxn,
  value,
  label,
  minWidth,
  integer
}: NumberPieceInputs): React.ReactElement {

  function onNewNum(event: any) {
    const parser = integer ? parseIntBetter : Number;
    const parsed = parser(event.target.value);
    if (!isNaN(parsed)) {
      changeFxn(parsed, name);
    } else {
      console.log(`Number update not used as it did not parse properly: ${event.target.value}`)
    }
  }

  return (
    <div>
      <TextField
        key={`${name}-number`}
        value={from_num(value, !integer)}
        label={label}
        onChange={onNewNum}
        size='small'
        style={{minWidth: minWidth || 200}}
      />
    </div>
  );
};

export function NumberPiece(
  key: NumberPieceInputs['name'],
  changeFxn: NumberPieceInputs['changeFxn'],
  value: NumberPieceInputs['value'] = null,
  label: NumberPieceInputs['label'] = '',
  minWidth: NumberPieceInputs['minWidth'] = 200,
  integer: NumberPieceInputs['integer'] = false
) {
  return <NumberPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    minWidth={minWidth}
    integer={integer}
  />
};

export function FloatPiece(
  key: NumberPieceInputs['name'],
  changeFxn: NumberPieceInputs['changeFxn'],
  value: NumberPieceInputs['value'] = null,
  label: NumberPieceInputs['label'] = '',
  minWidth: NumberPieceInputs['minWidth'] = 200,
) {
  return <NumberPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    minWidth={minWidth}
    integer={false}
  />
}

export function IntegerPiece(
  key: NumberPieceInputs['name'],
  changeFxn: NumberPieceInputs['changeFxn'],
  value: NumberPieceInputs['value'] = null,
  label: NumberPieceInputs['label'] = '',
  minWidth: NumberPieceInputs['minWidth'] = 200,
) {
  return <NumberPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    minWidth={minWidth}
    integer={true}
  />
}

interface SliderPieceInputs extends PieceBaseInputs<number> {
  min: number,
  max: number,
  stepSize?: number
}

export function SliderPieceRct({
  name, changeFxn, value, label,
  min,
  max,
  stepSize
}: SliderPieceInputs): React.ReactElement {
  return(
    <div key={name}>
      <InputLabel htmlFor={'slider-'+name} shrink>{label}</InputLabel>
      <Slider
        key={'slider-'+name}
        value={value}
        onChange={(event, newValue) => changeFxn(newValue as number, name)}
        min={min}
        max={max}
        step={stepSize}
        valueLabelDisplay="auto"
      />
    </div>
  );
}

export function SliderPiece(
  key: SliderPieceInputs['name'],
  changeFxn: SliderPieceInputs['changeFxn'],
  value: SliderPieceInputs['value'],
  label?: SliderPieceInputs['label'],
  min: SliderPieceInputs['min'] = 0.1,
  max: SliderPieceInputs['max'] = 20,
  stepSize: SliderPieceInputs['stepSize'] = undefined
) {
  return <SliderPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    min={min}
    max={max}
    stepSize={stepSize}
  />
}

export type NumericConstraint = ['exactly' | 'above', number | null, 'exactly' | 'below', number | null];
export const emptyNumericConstraintDef: NumericConstraint = ['exactly', 0, 'below', 0];

interface RangePieceInputs extends PieceBaseInputs<NumericConstraint> {

};

export function RangePieceRct({
  name,
  changeFxn,
  value,
  label = ''
}: RangePieceInputs): React.ReactElement {

  function updateSlot(newValue: 'exactly' | 'above' | 'below' | number | null, slot: number, current_full: NumericConstraint = value) {
    let next_full: NumericConstraint = [...current_full];
    next_full[slot] = newValue;
    changeFxn(next_full, name);
  };

  return(
    <div key={`${name}-range`}>
      {!!label ? <InputLabel htmlFor={`${name}-range-input`} shrink>{label}</InputLabel>: null}
      <Grid container direction='column' key={`${name}-range-input`} spacing={1} style={!!label ? {paddingLeft: 10} : {}}>
        <Grid item container spacing={1}>
          <Grid item xs={4}>
            <DropdownPieceRct
              name={name+'_lower_bound_type'}
              changeFxn={(newValue: string | null, k?: string) => {
                if (!!newValue && ['exactly','above'].includes(newValue)) {
                  updateSlot(newValue as 'exactly' | 'above', 0)
                }
              }}
              value={value[0]}
              label='From'
              options_in={['exactly','above']}
              minWidth={120}
            />
          </Grid>
          <Grid item xs={8}>
            <NumberPieceRct
              name={name+'_lower_value'}
              changeFxn={(newValue: number | null, k?: string) => {updateSlot(newValue, 1)}}
              value={value[1]}
              label='Min-value'
              minWidth={120}
              integer={false}
            />
          </Grid>
        </Grid>
        <Grid item container spacing={1}>
          <Grid item xs={4}>
            <DropdownPieceRct
              name={name+'_upper_bound_type'}
              changeFxn={(newValue: string | null, k?: string) => {
                if (!!newValue && ['exactly','below'].includes(newValue)) {
                  updateSlot(newValue as 'exactly' | 'below', 2)
                }
              }}
              value={value[2]}
              label='To'
              options_in={['exactly','below']}
              minWidth={120}
            />
          </Grid>
          <Grid item xs={8}>
            <NumberPieceRct
              name={name+'_upper_value'}
              changeFxn={(newValue: number | null, k?: string) => {updateSlot(newValue, 3)}}
              value={value[3]}
              label='Max-value'
              minWidth={120}
              integer={false}
            />
          </Grid>
        </Grid>
      </Grid>
    </div>
  );
}

export function rangePiece(
  key: RangePieceInputs['name'],
  changeFxn: RangePieceInputs['changeFxn'],
  value: RangePieceInputs['value'],
  label?: RangePieceInputs['label']
) {
  return <RangePieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
  />
}
