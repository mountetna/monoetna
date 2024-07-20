import React from 'react';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import {makeStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';

const useStyles = makeStyles((theme) => ({
  select: {
    alignItems: 'center',
    '&&& .MuiSelect-select': {
      paddingRight: '5px'
    }
  }
}));

function id(label: string) {
  return `${label}-${Math.random()}`;
}

const Selector = React.memo(
  ({
    canEdit,
    onSelect,
    label,
    name,
    placeholder,
    choiceSet,
    color='purple'
  }: {
    canEdit: boolean;
    onSelect: (name: string) => void;
    label: string;
    name: string;
    placeholder: string;
    choiceSet: string[]|string[][];
  }) => {
    const classes = useStyles();

    if (!canEdit) return <Typography>{name}</Typography>;
    const theme = useTheme();
    return (
      <FormControl variant="standard" className={classes.select}
        sx={{
          '&&& .MuiSelect-select': {
            color: color in theme.palette ? theme.palette[color].main : color
          }
        }}
      >
        <Select
          labelId={id(label)}
          variant='standard'
          autoWidth
          color={ color }
          IconComponent={null}
          disableUnderline={true}
          value={name}
          displayEmpty
          onChange={(e) => onSelect(e.target.value as string)}
        >
          <MenuItem disabled value=''><Typography color='gray'>{placeholder}</Typography></MenuItem>
          {choiceSet.sort().map((option: string|string[], index: number) => {
            let [ txt, value ] = [ option, option ];
            if (Array.isArray(option)) [ txt, value ] = option;
            return <MenuItem key={index} value={value}>
              {txt}
            </MenuItem>
          })}
        </Select>
      </FormControl>
    );
  }
);

export default Selector;
