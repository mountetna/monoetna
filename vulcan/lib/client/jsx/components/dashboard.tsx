import React, {useState} from 'react';
import 'regenerator-runtime/runtime';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

import {Workflow} from '../api_types';
import WorkflowsTable from './dashboard/workflows_table';
import FiguresTable from './dashboard/figures_table';

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

export default function Dashboard({project_name}: {project_name: string}) {
  const [selectedWorkflow, setSelectedWorkflow] = useState<Workflow | null>(
    null
  );
  const classes = useStyles();

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
          <WorkflowsTable
            project_name={project_name}
            onSelectWorkflow={(workflow) => setSelectedWorkflow(workflow)}
          />
        </Grid>
        <Grid item>
          <Typography variant='h6' className={classes.tableHeader}>
            Saved Figures
          </Typography>
        </Grid>
        <Grid item className={classes.workflows}>
          <FiguresTable
            project_name={project_name}
            workflowName={selectedWorkflow?.name}
          />
        </Grid>
      </Grid>
    </main>
  );
}
