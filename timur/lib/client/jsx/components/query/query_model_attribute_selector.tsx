import React, {useCallback, useState, useEffect, useMemo} from 'react';
import ModelMapGraphic from '../model_map/model_map_graphic';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';

import List from '@material-ui/core/List';
import ListItemButton from '@material-ui/core/ListItemButton';
import ListItemText from '@material-ui/core/ListItemText';
import ListItem from '@material-ui/core/ListItem';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectModelNames, selectTemplate} from 'etna-js/selectors/magma';
import MapSelector from './map_selector';

const useStyles = makeStyles((theme) => ({
  attributes: {
    maxHeight: '580px',
    overflow: 'scroll'
  },
  selection: {
    cursor: 'pointer',
    '&:hover': {
      textDecoration: 'underline'
    }
  }
}));

const QueryModelAttributeSelector = ({setModel, modelNames, modelName, setAttribute, attributeName}) => {
  const [ internalModelName, setInternalModel ] = useState(modelName);

  const updateAttribute = (modelName, attributeName) => {
    setAttribute(modelName, attributeName);
    setOpen(false);
  };
  const [ open, setOpen ] = useState(false);

  const readOnly = !setModel && !setAttribute;

  const openDialog = () => (!readOnly && setOpen(true));

  const classes = useStyles();

  return <>
    <Grid className={ !readOnly ? classes.selection : null } onClick={ openDialog }>
    {
      modelName
        ? <Typography component='span' color={ setModel ? 'secondary' : '#444' }>{modelName}</Typography>
        : <Typography component='span' color='red'>model_name</Typography>
    }
    <Typography component='span'>.</Typography>
    {
      attributeName
        ? <Typography component='span' color={ setAttribute ? 'secondary' : '#444' }>{attributeName}</Typography>
        : <Typography component='span' color='red'>attribute_name</Typography>
    }
    </Grid>
    <MapSelector
      open={open}
      onClose={() => setOpen(false)}
      setModel={setInternalModel}
      setAttribute={updateAttribute}
      modelNames={modelNames}
      attributeName={attributeName}
      modelName={internalModelName}/>
  </>
}

export default QueryModelAttributeSelector;
