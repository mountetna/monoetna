import React from 'react';
import Checkbox from "@material-ui/core/Checkbox";
import FormControlLabel, { FormControlLabelProps } from "@material-ui/core/FormControlLabel";
import { PieceBaseInputs } from './user_input_pieces';

interface CheckboxPieceInputs extends PieceBaseInputs<boolean> {
  disabled?: boolean
  labelPlacement?: FormControlLabelProps['labelPlacement']
}

export function CheckboxPiece(
  key: CheckboxPieceInputs['name'],
  changeFxn: CheckboxPieceInputs['changeFxn'],
  value: CheckboxPieceInputs['value'],
  label: CheckboxPieceInputs['label'] = '',
  disabled: CheckboxPieceInputs['disabled'] = false,
  labelPlacement: CheckboxPieceInputs['labelPlacement'] ='end') {

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
  labelPlacement = 'end'
}: CheckboxPieceInputs): React.ReactElement {

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