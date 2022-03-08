import React, {useState} from 'react';
import 'regenerator-runtime/runtime';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

import {Workflow} from '../api_types';
import WorkflowsCarousel from './dashboard/workflows_carousel';
import FiguresGrid from './dashboard/figures_grid';
import WorkflowControls from './dashboard/workflow_controls';

const useStyles = makeStyles((theme) => ({
  title: {
    padding: '10px 15px 5px',
    color: '#444'
  },
  workflows: {},
  tableHeader: {
    paddingLeft: '1rem',
    fontSize: '1.25rem',
    fontWeight: 'bold'
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
        <Grid item container justifyContent='space-between'>
          <Typography variant='h6' className={classes.tableHeader}>
            Available Workflows
          </Typography>
          <WorkflowControls
            workflow={selectedWorkflow}
            project_name={project_name}
          />
        </Grid>
        <Grid item className={classes.workflows}>
          <WorkflowsCarousel
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
          <FiguresGrid
            project_name={project_name}
            workflowName={selectedWorkflow?.name}
          />
        </Grid>
      </Grid>
    </main>
  );
}
