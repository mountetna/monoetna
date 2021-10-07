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

export default function Dashboard({project_name}) {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const classes = useStyles();

  const visitWorkflow = useCallback((workflow) => {
    invoke(pushLocation(`/${project_name}/workflow/${workflowName(workflow)}`));
  }, [invoke])

  const projectWorkflows = workflows ? workflows.filter( ({projects}) => projects.includes(project_name) ) : [];

  return (
    <main className='vulcan-dashboard'>
      <Grid container direction='column'>
        <Grid item container className={classes.title}><Typography variant='h5'>{project_name} workflows</Typography></Grid>
        <Grid item container className={classes.workflows}>
          {
            projectWorkflows.length ? projectWorkflows.map((w, ind) =>
                <Grid item key={ind}><Card
                  workflow={w}
                  onClick={() => visitWorkflow(w) }
              /></Grid>
                ) : <Grid item className={classes.none}><em>No workflows</em></Grid>
          }
        </Grid>
      </Grid>
    </main>
  );
}
