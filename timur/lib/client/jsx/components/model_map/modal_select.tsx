import React from 'react';

import TextField from '@material-ui/core/TextField';
import MenuItem from '@material-ui/core/MenuItem';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  popover: {
    zIndex: '30000 !important' as any // etna modal is 20000
  }
}));

export default function ModalSelect({
  onChange,
  options,
  id,
  label,
  value
}: {
  onChange: any;
  options: string[];
  id: string;
  label: string;
  value: string;
}) {
  const classes = useStyles();
  return (
    <TextField
      id={id}
      select
      value={value}
      label={label}
      SelectProps={{
        MenuProps: {
          PopoverClasses: {
            root: classes.popover
          }
        }
      }}
      onChange={(e: React.ChangeEvent<any>) => onChange(e.target.value)}
    >
      {options.sort().map((option: string, i: number) => (
        <MenuItem key={i} value={option}>
          {option}
        </MenuItem>
      ))}
    </TextField>
  );
}
