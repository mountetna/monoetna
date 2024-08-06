import React, {useState, useMemo, useContext, useCallback} from 'react';
import IconButton from '@mui/material/IconButton';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import ChevronRightIcon from '@mui/icons-material/ChevronRight';

import { makeStyles } from '@mui/styles';

const useStyles = makeStyles((theme) => ({
  chevron: {
    marginRight: '5px'
  }
}));

const QueryChevron = ({fold=true,setFold,disabled}) => {
  const classes = useStyles();
  return (
    <IconButton className={classes.chevron} disabled={disabled} size='small' onClick={ setFold ? (() => setFold(!fold)) : null } >
      { fold ? <ChevronRightIcon fontSize='small'/> : <ExpandMoreIcon fontSize='small'/> }
    </IconButton>
  );
};

export default QueryChevron;
