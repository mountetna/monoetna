import React, {useCallback, useContext, useState, useEffect} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';

import {VulcanContext} from '../contexts/vulcan_context';
import Card from '../components/dashboard/card';
import {workflowName} from "../selectors/workflow_selectors";
import SelectInput from 'etna-js/components/inputs/select_input'

import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import { json_get } from 'etna-js/utils/fetch';

const useStyles = makeStyles( theme => ({
  title: {
    padding: '10px 0px',
    color: '#444'
  },
  figures: {
    padding: '15px'
  },
  none: {
    color: '#f44'
  }
}));

export default function FigureList({project_name}) {
  const invoke = useActionInvoker();
  const [ figures, setFigures ] = useState(null);

  const classes = useStyles();

  useEffect( () => {
    json_get(
      `/api/${project_name}/figures`
    ).then(
      ({figures}) => setFigures(figures)
    )
  }, []);


  return (
    <main className={classes.figures}>
      <Grid container direction='column'>
        <Grid item container className={classes.title}><Typography variant='h5'>{project_name} figures</Typography></Grid>
        {
          figures ? figures.map(
            (figure,i) => <Figure figure={figure}/>
          ) : "No figures"
        }
        <Grid><Button onClick={showCreateFigure} variant="text">Create Figure</Button></Grid>
      </Grid>
    </main>
  );
}
