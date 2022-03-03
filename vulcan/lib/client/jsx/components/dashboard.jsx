import React, {useCallback, useContext, useState} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';

import {VulcanContext} from '../contexts/vulcan_context';
import WorkflowsTable from './dashboard/workflows_table';
import {workflowName} from '../selectors/workflow_selectors';
import FiguresTable from './dashboard/figures_table';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  title: {
    padding: '10px 15px 5px',
    color: '#444'
  },
  workflows: {
    padding: '15px'
  },
  tableHeader: {
    paddingLeft: '1rem',
    fontSize: '1rem',
    color: '#de5833'
  }
}));

export default function Dashboard({project_name}) {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const classes = useStyles();

  const visitWorkflow = useCallback(
    (workflow) => {
      invoke(
        pushLocation(`/${project_name}/workflow/${workflowName(workflow)}`)
      );
    },
    [invoke]
  );

  const projectWorkflows = workflows
    ? workflows.filter(({projects}) => projects.includes(project_name))
    : [];

  return (
    <main className='vulcan-dashboard'>
      <Grid container direction='column'>
        <Grid item container className={classes.title}>
          <Typography variant='h5'>{project_name}</Typography>
        </Grid>
        <Grid item>
          <Typography variant='h6' className={classes.tableHeader}>
            Available Workflows
          </Typography>
        </Grid>
        <Grid item className={classes.workflows}>
          <WorkflowsTable project_name={project_name} />
        </Grid>
        <Grid item>
          <Typography variant='h6' className={classes.tableHeader}>
            Saved Figures
          </Typography>
        </Grid>
        <Grid item className={classes.workflows}>
          <FiguresTable project_name={project_name} />
        </Grid>
      </Grid>
    </main>
  );
}
