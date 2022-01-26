import React, {useCallback, useContext, useState} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';

import {VulcanContext} from '../contexts/vulcan_context';
import Card from '../components/dashboard/card';
import {workflowName} from "../selectors/workflow_selectors";
import SelectInput from 'etna-js/components/inputs/select_input'

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles( theme => ({
  title: {
    padding: '10px 15px 5px',
    color: '#444'
  },
  workflows: {
    padding: '15px'
  },
  none: {
    color: '#f44'
  }
}));

export default function FigureList({project_name}) {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {figures} = state;

  const classes = useStyles();

  return (
    <main className='vulcan-figures'>
      <Grid container direction='column'>
        <Grid item container className={classes.title}><Typography variant='h5'>{project_name} figures</Typography></Grid>
      </Grid>
    </main>
  );
}
