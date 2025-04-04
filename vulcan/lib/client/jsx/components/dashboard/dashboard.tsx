import React, {useContext, useEffect, useMemo, useState} from 'react';
import 'regenerator-runtime/runtime';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

import {Workflow, WorkspacesResponse} from '../../api_types';
import WorkflowsCarousel from './workflows_carousel';
import WorkspacesGrid from './workspace_control/workspaces_grid';
import WorkspacesControls from './workspace_control/workspaces_controls';
import WorkflowControls from './workflow_control/workflow_controls';
import ProjectHeader from 'etna-js/components/project-header';
import { VulcanContext } from '../../contexts/vulcan_context';


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
  headerTitle: {}
}));

export default function Dashboard({project_name}: {project_name: string}) {
  const [selectedWorkflow, setSelectedWorkflow] = useState<Workflow | null>(
    null
  );
  const classes = useStyles();

  const [searchString, setSearchString] = useState('');
  const [tags, setTags] = useState<string[]>(['public']);

  let {state} = useContext(VulcanContext);
  const {workflows, workspaces} = useMemo(() => {
    return state
  }, [state.workflows, state.workspaces])

  function trim_workspaces_to_project(wflows, p_name, wpaces) {
    const projectWorkflowIds = wflows.filter(w => w.project_name == p_name).map(w => w.id)
    return wpaces.filter((w) => projectWorkflowIds.includes(w.workflow_id))
  }
  const projectWorkspaces = useMemo(() => trim_workspaces_to_project(workflows, project_name, workspaces),
  [workflows, workspaces, project_name]);

  return (
    <main className='vulcan-dashboard'>
      <Grid container direction='column'>
        <ProjectHeader project_name={ project_name } className={classes.title}/>
        <Grid item container className={classes.tableHeader}>
          <Grid item xs={2}>
            <Typography
              color='secondary'
              variant='h6'
              className={classes.headerTitle}
            >
              Workflows
            </Typography>
          </Grid>
          <Grid item xs={10}>
            <WorkflowControls
              workflow={selectedWorkflow}
              project_name={project_name}
              workspaces={workspaces}
            />
          </Grid>
        </Grid>
        <WorkflowsCarousel
          project_name={project_name}
          onSelectWorkflow={(workflow) => setSelectedWorkflow(workflow)}
        />
        <Grid item container className={classes.tableHeader}>
          <Grid item xs={2}>
            <Typography
              color='secondary'
              variant='h6'
              className={classes.headerTitle}
            >
              Workspaces
            </Typography>
          </Grid>
          <Grid item xs={10}>
            <WorkspacesControls
              setSearchString={setSearchString}
              setTags={setTags}
              workspaces={projectWorkspaces}
              tags={tags}
              searchString={searchString}
            />
          </Grid>
        </Grid>
        <WorkspacesGrid
          project_name={project_name}
          workflowId={selectedWorkflow?.id}
          tags={tags}
          searchString={searchString}
          allWorkflows={workflows}
          allWorkspaces={projectWorkspaces}
        />
      </Grid>
    </main>
  );
}
