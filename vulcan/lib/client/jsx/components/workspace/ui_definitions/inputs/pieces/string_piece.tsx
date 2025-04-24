import React from "react";
import TextField from "@material-ui/core/TextField";
import { PieceBaseInputs } from "./user_input_pieces";

interface StringPieceInputs extends PieceBaseInputs<string> {
  minWidth?: number;
  multiline?: boolean;
  fullWidth?: boolean;
}

export function StringPieceRct({
  name, changeFxn, value,
  label = '',
  minWidth = 150,
  multiline = undefined,
  fullWidth = undefined,
}: StringPieceInputs): React.ReactElement {
  return (
  <div key={name}>
    <TextField
      value={value}
      multiline={multiline}
      fullWidth={fullWidth}
      label={label}
      InputLabelProps={{ shrink: true }}
      onChange={(event) => changeFxn(event.target.value, name)}
      size="small"
      style={{minWidth: minWidth || 200}}
    />
  </div>
  );
};

export function StringPiece(
  key: StringPieceInputs['name'],
  changeFxn: StringPieceInputs['changeFxn'],
  value: StringPieceInputs['value'],
  label: StringPieceInputs['label'] = '',
  minWidth: StringPieceInputs['minWidth'] = 150,
  multiline: StringPieceInputs['multiline'] = undefined,
  fullWidth: StringPieceInputs['fullWidth'] = undefined,
) {
  return <StringPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    minWidth={minWidth}
    multiline={multiline}
    fullWidth={fullWidth}
  />
};