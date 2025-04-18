import React, {useState, useMemo, useContext, useCallback} from 'react';
import IconButton from '@material-ui/core/IconButton';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ChevronRightIcon from '@material-ui/icons/ChevronRight';

import { makeStyles } from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  chevron: {
    marginRight: '5px'
  }
}));

const QueryChevron = ({fold=true,setFold=(f: boolean) => {},disabled=false}:{
  fold?: boolean;
  setFold?: (f: boolean) => void;
  disabled?: boolean;
}) => {
  const classes = useStyles();
  return (
    <IconButton className={classes.chevron} disabled={disabled} size='small' onClick={ () => setFold(!fold) }>
      { fold ? <ChevronRightIcon fontSize='small'/> : <ExpandMoreIcon fontSize='small'/> }
    </IconButton>
  );
};

export default QueryChevron;
