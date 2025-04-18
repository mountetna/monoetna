import React from 'react';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import {makeStyles, useTheme} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import { PaletteColor } from '@material-ui/core/styles/createPalette';

const useStyles = (color: string) => makeStyles((theme) => ({
  select: {
    alignItems: 'center',
    '&&& .MuiSelect-select': {
      padding: '6px 0px',
      paddingRight: '5px',
      color: color in theme.palette ? (theme.palette[color as keyof typeof theme.palette] as PaletteColor).main : color
    }
  }
}));

function id(label: string) {
  return `${label}-${Math.random()}`;
}

const Selector = ({
  canEdit,
  onSelect,
  label='select',
  name,
  placeholder,
  choiceSet,
  color='purple'
}: {
  canEdit: boolean;
  onSelect: (name: string) => void;
  label?: string;
  name: string;
  placeholder?: string;
  choiceSet: (string|[string,boolean])[];
  color?: string;
}) => {
  const classes = useStyles(color)();

  if (!canEdit) return <Typography>{name}</Typography>;
  const theme = useTheme();
  return (
    <FormControl variant="standard" className={classes.select}>
      <Select
        labelId={id(label)}
        variant='standard'
        autoWidth
        IconComponent='span'
        disableUnderline={true}
        value={name}
        displayEmpty
        onChange={(e) => onSelect(e.target.value as string)}
      >
        { placeholder && <MenuItem disabled value=''><Typography style={{ color:'gray'}} >{placeholder}</Typography></MenuItem> }
        {choiceSet.sort().map((option: string|[string, boolean], index: number) => {
          let [ txt, value ]: [string,any] = Array.isArray(option) ? option : [ option, option ];
          return <MenuItem key={index} value={value}>
            {txt}
          </MenuItem>
        })}
      </Select>
    </FormControl>
  );
}

export default Selector;
