import React from 'react';
import Typography from '@mui/material/Typography';
import Tooltip from '@mui/material/Tooltip';
import IconButton from '@mui/material/IconButton';
import AddIcon from '@mui/icons-material/Add';
import RestartAltIcon from '@mui/icons-material/RestartAlt';

import { makeStyles } from '@mui/styles';

const useStyles = makeStyles((theme) => ({
  folded: {
    fontStyle: 'italic',
    paddingLeft: '10px',
    cursor: 'pointer'
  }
}));

const QueryClauseSummaryControls = ({fold, setFold, itemName, numItems, addHandler, removeHandler}) => {
  const classes = useStyles();
  return (<>
    {
      fold
        ? <Typography
            component='span'
            onClick={ () => setFold(!fold) }
            className={classes.folded}>
            { numItems || 'no' } {itemName}{ numItems == 1 ? '' : 's' }
          </Typography>
      : <>
          <Tooltip title={`Add ${itemName}`} aria-label={`Add ${itemName}`}>
            <IconButton
              size='small'
              onClick={ addHandler }
              color='primary'>
              <AddIcon fontSize='small'/>
            </IconButton>
          </Tooltip>
          <Tooltip title={`Remove all ${itemName}s`} aria-label={`Remove all ${itemName}s`}>
            <IconButton size='small' onClick={removeHandler} color='primary'>
              <RestartAltIcon fontSize='small'/>
            </IconButton>
          </Tooltip>
        </>
    }
  </>);
};

export default QueryClauseSummaryControls;
