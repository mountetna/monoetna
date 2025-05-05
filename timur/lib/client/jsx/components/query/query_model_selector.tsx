import React, {useState} from 'react';
import InputLabel from '@material-ui/core/InputLabel';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import {makeStyles} from '@material-ui/core/styles';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';

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

const QueryModelSelector = ({setModel, modelNames, modelName}:{
  setModel?: (modelName: string) => void;
  modelNames: string[];
  modelName: string;
}) => {
  const [ open, setOpen ] = useState(false);

  const updateModel = (newModelName: string) => {
    if (setModel) setModel(newModelName);
    setOpen(false);
  };

  const readOnly = !setModel;

  const openDialog = () => (!readOnly && setOpen(true));

  const classes = useStyles();

  return <>
    <Grid className={ !readOnly ? classes.selection : undefined } onClick={ openDialog }>
    {
      modelName
        ? <Typography component='span' color={ setModel ? 'secondary' : 'inherit' }>{modelName}</Typography>
        : <Typography component='span' color='error'>model_name</Typography>
    }
    </Grid>
    <MapSelector
      open={open}
      onClose={() => setOpen(false)}
      setModel={updateModel}
      modelNames={modelNames}
      modelName={modelName}/>
  </>
}

export default QueryModelSelector;
