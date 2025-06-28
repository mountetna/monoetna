import React, {useContext, useEffect, useMemo, useState} from 'react';
import 'regenerator-runtime/runtime';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

import {Workflow} from '../../api_types';
import WorkflowsCarousel from './workflows_carousel';
import WorkspacesGrid from './workspace_control/workspaces_grid';
import WorkspacesControls from './workspace_control/workspaces_controls';
import WorkflowControls from './workflow_control/workflow_controls';
import ProjectHeader from 'etna-js/components/project-header';
import { VulcanContext } from '../../contexts/vulcan_context';
import {useFeatureFlag} from "etna-js/hooks/useFeatureFlag";
import useUserHooks from '../../contexts/useUserHooks';


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

  if (!useFeatureFlag('vulcan')) {
    return <div style={{padding: '10px'}}>
      <Typography variant='h3'>
        Vulcan is under development
      </Typography>
      <Typography>
        Access is restricted to a small group of testers at this stage, and it appears that you are not one of these testers.
      </Typography>
      <Typography>
        Contact the engineering team if you believe you are seeing this message in error.
      </Typography>
    </div>
  }

  const [selectedWorkflow, setSelectedWorkflow] = useState<Workflow | null>(
    null
  );
  const classes = useStyles();

  const [searchString, setSearchString] = useState('');
  const [tags, setTags] = useState<string[]>(['highlighted']);

  let {state} = useContext(VulcanContext);
  const {workflows, workspaces} = useMemo(() => {
    return state
  }, [state.workflows, state.workspaces])

  const {canEdit} = useUserHooks();
  const visibleWorkspaces = useMemo(() => {
    const projectWorkflowIds = workflows.filter(w => w.project_name == project_name).map(w => w.id)
    const projectWorkspaces = workspaces.filter((w) => projectWorkflowIds.includes(w.workflow_id))
    return projectWorkspaces.filter((w) => canEdit(w) || (w.tags || []).includes('published'))
  }, [workflows, workspaces, project_name]);

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
              workspaces={visibleWorkspaces}
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
          accessibleWorkspaces={visibleWorkspaces}
        />
      </Grid>
    </main>
  );
}
