import React, { useState, useEffect, useCallback } from 'react';
import Grid from '@material-ui/core/Grid';

import {makeStyles} from '@material-ui/core/styles';
import { getDocuments } from 'etna-js/api/magma_api';

const useStyles = makeStyles( theme => ({
  form: {
    border: '1px solid #ccc'
  }
}));

const RedcapAttribute = ({name, value}) => {
  return <Grid item container>
    <Grid item>{JSON.stringify(value)}</Grid>
  </Grid>
}

const RedcapScript = ({script}) => {
  const { attributes, each } = script;

  return <Grid container>
    {
      Object.keys(attributes).map(
        att_name => <Grid key={att_name} item>
          <RedcapAttribute name={att_name} value={attributes[att_name]}/>
        </Grid>
      )
    }
  </Grid>
}

const RedcapModel = ({config,modelName}) => {
  const classes = useStyles();
  const { each, invert, scripts } = config;
  return <Grid container>
    <Grid item xs={2}>{ modelName }</Grid>
    <Grid item xs={2}>{ JSON.stringify(each) }</Grid>
    <Grid item xs={1}>{ invert }</Grid>
    <Grid item xs={7}>{ scripts.map(
      (script,i) => <RedcapScript key={i} script={script}/>
    )}</Grid>
  </Grid>
};

const RedcapForm = ({config}:{
  config: any
}) => {
  const classes = useStyles();

  const modelNames = Object.keys(config);

  useEffect( () => {
    getDocuments({
      model_name: 'all',
      record_names: [],
      attribute_names: 'all',
    }, fetch).then( response => console.log({response}) ).catch( e => e.then( r => console.log({r}) ) )
  });

  return <Grid container direction='column' className={classes.form}>
    <Grid item> model_name </Grid>
    {
      modelNames.map( modelName => <RedcapModel key={modelName} modelName={modelName}
        config={config[modelName]}
      />)
    }
  </Grid>
}

export default RedcapForm;
