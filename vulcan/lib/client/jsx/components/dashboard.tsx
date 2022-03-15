import React, {useState} from 'react';
import 'regenerator-runtime/runtime';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

import {Workflow} from '../api_types';
import WorkflowsCarousel from './dashboard/workflows_carousel';
import FiguresGrid from './dashboard/figures_grid';
import FiguresControls from './dashboard/figures_controls';
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
    borderBottom: '1px solid #eee',
    borderTop: '1px solid #eee',
    alignItems: 'center',
    minHeight: '88px'
  },
  headerTitle: {
    flex: '0 0 330px'
  }
}));

export default function Dashboard({project_name}: {project_name: string}) {
  const [selectedWorkflow, setSelectedWorkflow] = useState<Workflow | null>(
    null
  );
  const classes = useStyles();

  const [searchString, setSearchString] = useState('');
  const [tags, setTags] = useState<string[]>(['public']);

  return (
    <main className='vulcan-dashboard'>
      <Grid container direction='column'>
        <Grid item container className={classes.title}>
          <Typography color='primary' variant='h5'>{project_name}</Typography>
        </Grid>
        <Grid item container className={classes.tableHeader}>
          <Typography color='secondary' variant='h6' className={classes.headerTitle}>
            Workflows
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
        <Grid item container className={classes.tableHeader}>
          <Typography color='secondary' variant='h6' className={classes.headerTitle}>
            Figures
          </Typography>
          <FiguresControls
            setSearchString={setSearchString}
            setTags={setTags}
            project_name={project_name} />
        </Grid>
        <Grid item className={classes.workflows}>
          <FiguresGrid
            project_name={project_name}
            workflowName={selectedWorkflow?.name}
            tags={tags}
            searchString={searchString}
          />
        </Grid>
      </Grid>
    </main>
  );
}
