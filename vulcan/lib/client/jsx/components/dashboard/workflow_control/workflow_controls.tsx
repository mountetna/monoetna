import React from 'react';
import Grid from '@material-ui/core/Grid';
import WorkflowCreateButtonModal from './workflow_creation';
import WorkspaceCreateButtonModal from './workspace_creation';
import { VulcanState } from '../../../reducers/vulcan_reducer';
import useUserHooks from '../../../contexts/useUserHooks';

export default function WorkflowControls({
  project_name,
  workflow,
  workspaces,
  workflows,
}: {
  project_name: string;
  workflow: VulcanState['workflow'] | null;
  workspaces: VulcanState['workspaces'];
  workflows: VulcanState['workflows']
}) {

  const {superuser} = useUserHooks();

  return (
    <Grid
      container
      alignItems='center'
    >
      <Grid item>
        <WorkspaceCreateButtonModal projectName={project_name} workflow={workflow} workspaces={workspaces} />
      </Grid>
      {superuser &&
      <Grid item>
        <WorkflowCreateButtonModal projectName={project_name} workflows={workflows}/>
      </Grid>}
    </Grid>
  );
}
