import React, {useCallback, useState, useEffect, useMemo} from 'react';
import ModelMapGraphic from '../model_map/model_map_graphic';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';
import { makeStyles } from '@material-ui/core/styles';

import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import Typography from '@material-ui/core/Typography';

import {useReduxState} from 'etna-js/hooks/useReduxState';
import {Attribute} from 'etna-js/models/magma-model';
import {selectModelNames, selectTemplate} from 'etna-js/selectors/magma';
import ModelAttributesTable, {toggleSelection} from '../model_map/model_attributes_table';

const useStyles = makeStyles((theme) => ({
  attributes: {
    maxHeight: '600px',
    width: '480px',
    overflow: 'hidden',
    flexWrap: 'nowrap',
    padding: '0 25px 25px 0'
  },
  selection: {
    cursor: 'pointer',
    '&:hover': {
      textDecoration: 'underline'
    }
  }
}));

const MapSelector = ({setModel, open, onClose, modelNames=[], modelName, setAttribute, setAttributes, filterAttributes, attributeName}:{
  setAttribute?: (modelName: string, attributeName: string) => void;
  setAttributes?: (attributeNames: string[]) => void;
  setModel?: (modelName: string) => void;
  modelName: string;
  modelNames?: string[];
  attributeName?: string;
  open: boolean;
  filterAttributes?: (attribute: Attribute) => boolean;
  onClose: () => void;
}) => {
  const [ width, height ] = [ 600, 600 ];

  let {model_names, templates} = useReduxState((state: any) => {
    let model_names: string[] = selectModelNames(state);
    return {
      model_names,
      templates: Object.fromEntries(
        model_names.map((model_name) =>
          [ model_name, selectTemplate(state, model_name) ]
        )
      )
    };
  });

  const updateModel = (modelName: string) => {
    if (setModel) setModel(modelName);
    if (setAttributes) setSelected({});
  };

  let disabled = model_names.filter( (model_name: string) => modelNames.length && !modelNames.includes(model_name) );

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

  const handleClose = useCallback(
    () => {
      setSelected({});
      onClose();
    }, [ onClose ]
  );

  const handleSetAttributes = useCallback(
    () => {
      if (setAttributes) setAttributes(Object.keys(selected));
      handleClose();
    }, [ selected, handleClose ]
  );

  const setAttributeHandler = useCallback(
    attName => {
      if (setAttribute) setAttribute(modelName, attName);
      if (setAttributes) setSelected(toggleSelection(selected, attName));
    }, [ modelName, selected ]
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
            setAttribute={ setAttributeHandler }
            selected={ setAttributes ? selected : {} }
            setSelected={ setAttributes ? setSelected : undefined }
          />
        }
      </Grid>
  </Dialog>
}

export default MapSelector;
