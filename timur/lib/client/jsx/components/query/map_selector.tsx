import React, {useCallback, useState, useEffect, useMemo} from 'react';
import ModelMapGraphic from '../model_map/model_map_graphic';
import Dialog from '@mui/material/Dialog';
import DialogTitle from '@mui/material/DialogTitle';
import DialogContent from '@mui/material/DialogContent';
import { makeStyles } from '@mui/styles';

import Grid from '@mui/material/Grid';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';

import List from '@mui/material/List';
import ListItemButton from '@mui/material/ListItemButton';
import ListItemText from '@mui/material/ListItemText';
import ListItem from '@mui/material/ListItem';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectModelNames, selectTemplate} from 'etna-js/selectors/magma';
import ModelAttributesTable from '../model_map/model_attributes_table';

const useStyles = makeStyles((theme) => ({
  attributes: {
    maxHeight: '600px',
    width: '480px',
    overflow: 'hidden',
    padding: '0 25px 25px 0'
  },
  selection: {
    cursor: 'pointer',
    '&:hover': {
      textDecoration: 'underline'
    }
  }
}));

const MapSelector = ({setModel, open, onClose, modelNames=[], modelName, setAttribute, setAttributes, filterAttributes, attributeName}) => {
  const [ width, height ] = [ 600, 600 ];

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

  const updateModel = modelName => {
    if (setModel) setModel(modelName);
    if (setAttributes) setSelected({});
  };

  let disabled = model_names.filter( model_name => modelNames.length && !modelNames.includes(model_name) );

  const classes = useStyles();

  const filteredTemplate = useMemo( () => {
    if (!templates || !templates[modelName]) return null;

    const currentTemplate = templates[modelName];
    return {
      ...currentTemplate,
      attributes: Object.fromEntries(
        Object.keys( currentTemplate.attributes).filter(
          att_name => filterAttributes ? filterAttributes(currentTemplate.attributes[att_name]) : true
        ).map( att_name => [ att_name, currentTemplate.attributes[ att_name ] ])
      )
    };
  }, [ modelName, templates ] );

  const canSetAttribute = (setAttribute || setAttributes) && modelName && filteredTemplate

  const [ selected, setSelected ] = useState({});

  const handleSetAttributes = useCallback(
    () => {
      setAttributes(Object.keys(selected));
      handleClose();
    }, [ selected, handleClose ]
  );

  const handleClose = useCallback(
    () => {
      setSelected({});
      onClose();
    }, [ onClose ]
  );

  const attributesSelected = Object.keys(selected).length > 0;

  return <Dialog open={open} onClose={handleClose} maxWidth='lg'>
    <DialogTitle><Grid container justifyContent='space-between'>
        <span>Select model { setAttribute ? 'and attribute' : setAttributes ? 'and attributes' : ''}</span>
        { setAttributes && <Button onClick={ handleSetAttributes } disabled={ !attributesSelected }>Select Attributes</Button> }
    </Grid></DialogTitle>
      <Grid container>
        <ModelMapGraphic
          width={width}
          height={height}
          selected_models={[ modelName ].filter(_ => _)}
          disabled_models={disabled}
          handler={updateModel}
        />
        {
          canSetAttribute && <ModelAttributesTable
            className={classes.attributes}
            template={filteredTemplate}
            showHiddenAttributes={false}
            columns={ { type: true, attribute: true } }
            setAttribute={ setAttribute ? attName => setAttribute(modelName, attName) : null }
            selected={ setAttributes ? selected : {} }
            setSelected={ setAttributes ? setSelected : null }
          />
        }
      </Grid>
  </Dialog>
}

export default MapSelector;
