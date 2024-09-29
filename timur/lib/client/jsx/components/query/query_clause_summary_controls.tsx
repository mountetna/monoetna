import React from 'react';
import Typography from '@material-ui/core/Typography';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';
import AddIcon from '@material-ui/icons/Add';
import RemoveCircleIcon from '@material-ui/icons/RemoveCircle';

import { makeStyles } from '@material-ui/core/styles';

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
              <RemoveCircleIcon fontSize='small'/>
            </IconButton>
          </Tooltip>
        </>
    }
  </>);
};

export default QueryClauseSummaryControls;
