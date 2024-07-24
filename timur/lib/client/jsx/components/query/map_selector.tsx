import React, {useCallback, useState, useEffect, useMemo} from 'react';
import ModelMapGraphic from '../model_map/model_map_graphic';
import Dialog from '@mui/material/Dialog';
import DialogTitle from '@mui/material/DialogTitle';
import DialogContent from '@mui/material/DialogContent';
import { makeStyles } from '@mui/styles';

import Grid from '@mui/material/Grid';
import Typography from '@mui/material/Typography';

import List from '@mui/material/List';
import ListItemButton from '@mui/material/ListItemButton';
import ListItemText from '@mui/material/ListItemText';
import ListItem from '@mui/material/ListItem';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectModelNames, selectTemplate} from 'etna-js/selectors/magma';

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

const MapSelector = ({setModel, options, modelName, setAttribute, attributeName}) => {
  const [ width, height ] = [ 600, 600 ];
  const [ internalModelName, setInternalModel ] = useState(modelName);
  console.log({internalModelName, modelName, attributeName});
  const updateModel = (modelName) => {
    if (setAttribute) {
      setInternalModel(modelName);
      return;
    }
    setModel(modelName);
    setOpen(false);
  };
  const updateAttribute = (modelName, attributeName) => {
    setAttribute(modelName, attributeName);
    setOpen(false);
  };
  const [ open, setOpen ] = useState(false);
  let {model_names, templates} = useReduxState((state) => {
    let model_names = selectModelNames(state);
    return {
      model_names,
      templates: Object.fromEntries(
        model_names.map((model_name) =>
          [ model_name, selectTemplate(state, model_name) ]
        )
      )
    };
  });

  let disabled = model_names.filter( model_name => options.length && !options.includes(model_name) );

  const readOnly = !setModel && !setAttribute;

  const openDialog = () => (!readOnly && setOpen(true));

  const classes = useStyles();

  return <Grid>
    <Grid className={ !readOnly ? classes.selection : null } onClick={ openDialog }>
    {
      modelName
        ? <Typography component='span' color={ setModel ? 'secondary' : '#444' }>{modelName}</Typography>
        : <Typography component='span' color='red'>model_name</Typography>
    }{
      attributeName && <>
      <Typography component='span'>.</Typography>
      { 
        attributeName
          ? <Typography component='span' color={ setAttribute ? 'secondary' : '#444' }>{attributeName}</Typography>
          : <Typography component='span' color='red'>attribute_name</Typography>
      }
      </>
    }
    </Grid>
    <Dialog open={open} onClose={() => setOpen(false)} maxWidth='lg'>
      <DialogTitle>Select model { setAttribute ? 'and attribute' : ''}</DialogTitle>
        <Grid container>
          <ModelMapGraphic
            width={width}
            height={height}
            selected_models={[ internalModelName || modelName ].filter(_ => _)}
            disabled_models={disabled}
            handler={updateModel}
          />
          <Grid item maxHeight={600}>
          {
            internalModelName && templates[ internalModelName ] &&
            <List className={classes.attributes}>
              {
                Object.keys(templates[ internalModelName ].attributes).map(
                 attName =>
                   <ListItem disablePadding>
                     <ListItemButton onClick={ () => updateAttribute(internalModelName, attName) } selected={ attributeName == attName }>
                       <ListItemText primary={ attName }/>
                     </ListItemButton>
                   </ListItem>
                 )
               }
            </List>
          }
          </Grid>
        </Grid>
    </Dialog>
  </Grid>
}

export default MapSelector;
