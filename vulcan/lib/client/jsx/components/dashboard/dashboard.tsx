import React, {useContext, useEffect, useMemo, useState} from 'react';
import 'regenerator-runtime/runtime';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

import {Workflow} from '../../api_types';
import WorkspacesGrid from './workspaces_grid';
import WorkspacesFilterBar from './workspaces_filter_bar';
import ProjectHeader from 'etna-js/components/project-header';
import { VulcanContext } from '../../contexts/vulcan_context';
import {useFeatureFlag} from "etna-js/hooks/useFeatureFlag";
import {getProjects} from "etna-js/api/janus-api";
import useUserHooks from '../../contexts/useUserHooks';
import WorkflowCreateButtonModal from './workflow_creation';
import WorkspaceCreateButtonModal from './workspace_creation';


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

  const [privateSearchString, setPrivateSearchString] = useState('');
  const [privateTags, setPrivateTags] = useState<string[]>([]);
  const [publishedSearchString, setPublishedSearchString] = useState('');
  const [publishedTags, setPublishedTags] = useState<string[]>([]);

  // ToDo: Move this to VulcanContext or the etna-js ProjectHeader component itself.
  const [projectNameShow, setProjectNameShow] = useState(project_name);
  useEffect(() => {
    getProjects().then((r: any) => {
      const project_name_full: string[] = r.filter((p: any) => p.project_name==project_name).map((p: any) => p.project_name_full)
      if (project_name_full.length > 0) {
        setProjectNameShow(project_name_full[0])
      }
    })
  }, [project_name])

  let {state} = useContext(VulcanContext);
  const {workflows, workspaces} = useMemo(() => {
    return state
  }, [state.workflows, state.workspaces])

  const {canEdit, superuser} = useUserHooks();
  const visibleWorkspaces = useMemo(() => {
    const projectWorkflowIds = workflows.filter(w => w.project_name == project_name).map(w => w.id)
    const projectWorkspaces = workspaces.filter((w) => projectWorkflowIds.includes(w.workflow_id))
    return projectWorkspaces.filter((w) => canEdit(w) || (w.tags || []).includes('published'))
  }, [workflows, workspaces, project_name]);

  const personalWorkspaces = useMemo(() => {
    return visibleWorkspaces.filter((w) => canEdit(w))
  }, [workflows, workspaces, project_name]);

  const publishedWorkspaces = useMemo(() => {
    return visibleWorkspaces.filter((w) => (w.tags || []).includes('published'))
  }, [workflows, workspaces, project_name]);

  return (
    <main className='vulcan-dashboard'>
      <Grid container direction='column'>
        <ProjectHeader project_name={ project_name } project_name_full={projectNameShow}  className={classes.title}/>
        <Grid item container className={classes.tableHeader}>
          <Grid item xs={2}>
            <Typography
              color='secondary'
              variant='h6'
              className={classes.headerTitle}
            >
              Your Workspaces
            </Typography>
          </Grid>
          <Grid item xs={2}>
            <WorkspaceCreateButtonModal
              projectName={project_name}
              workspaces={workspaces}
            />
          </Grid>
          <Grid item xs={8}>
            <WorkspacesFilterBar
              setSearchString={setPrivateSearchString}
              setTags={setPrivateTags}
              workspaces={visibleWorkspaces}
              tags={privateTags}
              searchString={privateSearchString}
            />
          </Grid>
        </Grid>
        <WorkspacesGrid
          project_name={project_name}
          workflowId={selectedWorkflow?.id}
          tags={privateTags}
          searchString={privateSearchString}
          allWorkflows={workflows}
          accessibleWorkspaces={personalWorkspaces}
        />
        <Grid item container className={classes.tableHeader}>
          <Grid item xs={2}>
            <Typography
              color='secondary'
              variant='h6'
              className={classes.headerTitle}
            >
              Published Workspaces
            </Typography>
          </Grid>
          <Grid item xs={10}>
            <WorkspacesFilterBar
              setSearchString={setPublishedSearchString}
              setTags={setPublishedTags}
              workspaces={publishedWorkspaces}
              tags={publishedTags}
              searchString={publishedSearchString}
            />
          </Grid>
        </Grid>
        <WorkspacesGrid
          project_name={project_name}
          workflowId={selectedWorkflow?.id}
          tags={publishedTags}
          searchString={publishedSearchString}
          allWorkflows={workflows}
          accessibleWorkspaces={publishedWorkspaces}
        />
        {superuser &&
        <Grid item container className={classes.tableHeader}>
          <Grid item xs={2}>
            <Typography
              color='secondary'
              variant='h6'
              className={classes.headerTitle}
            >
              Manage Workflows
            </Typography>
          </Grid>
          <Grid item xs={10}>
            <WorkflowCreateButtonModal projectName={project_name} workflows={workflows}/>
          </Grid>
        </Grid>}
      </Grid>
    </main>
  );
}
