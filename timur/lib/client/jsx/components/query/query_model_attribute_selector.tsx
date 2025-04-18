import React, {useCallback, useState, useEffect, useMemo} from 'react';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';

import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectModelNames, selectTemplate} from 'etna-js/selectors/magma';
import {makeStyles, useTheme} from '@material-ui/core/styles';
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

const QueryModelAttributeSelector = ({setModel, modelNames, modelName, setAttribute, attributeName}:{
  modelName: string;
  attributeName?: string;
  setModel?: (modelName: string) => void;
  modelNames: string[];
  setAttribute?: (modelName: string, attributeName: string) => void;
}) => {
  const [ internalModelName, setInternalModel ] = useState(modelName);

  const updateAttribute = (modelName: string, attributeName: string) => {
    if (setAttribute) setAttribute(modelName, attributeName);
    setOpen(false);
  };
  const [ open, setOpen ] = useState(false);

  const readOnly = !setModel && !setAttribute;

  const openDialog = () => (!readOnly && setOpen(true));

  const classes = useStyles();

  const theme = useTheme();

  return <>
    <Grid className={ !readOnly ? classes.selection : undefined } onClick={ openDialog }>
    {
      modelName
        ? <Typography component='span' style={{ color: setModel ? theme.palette.secondary.main : '#444' }}>{modelName}</Typography>
        : <Typography component='span' color='error'>model_name</Typography>
    }
    <Typography component='span'>.</Typography>
    {
      attributeName
        ? <Typography component='span' style={{ color: setAttribute ? theme.palette.secondary.main : '#444' }}>{attributeName}</Typography>
        : <Typography component='span' color='error'>attribute_name</Typography>
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
