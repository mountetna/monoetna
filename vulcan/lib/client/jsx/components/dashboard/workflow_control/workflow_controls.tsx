import React, {useCallback} from 'react';

import Tooltip from '@material-ui/core/Tooltip';
import Button from '@material-ui/core/Button';
import InsertChartIcon from '@material-ui/icons/InsertChart';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';

import {pushLocation} from 'etna-js/actions/location_actions';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import WorkflowCreateButtonModal from './workflow_creation';
import WorkspaceCreateButtonModal from './workspace_creation';
import { VulcanState } from '../../../reducers/vulcan_reducer';

const useStyles = makeStyles((theme) => ({
  controls: {
  }
}));

export default function WorkflowControls({
  project_name,
  workflow,
  workspaces,
}: {
  project_name: string;
  workflow: VulcanState['workflow'];
  workspaces: VulcanState['workspaces'];
}) {
  const classes = useStyles();
  const invoke = useActionInvoker();

  return (
    <Grid
      container
      spacing={2}
    >
      <Grid item>
        <WorkflowCreateButtonModal projectName={project_name}/>
      </Grid>
      <Grid item>
        <WorkspaceCreateButtonModal projectName={project_name} workflow={workflow} workspaces={workspaces} />
      </Grid>
    </Grid>
  );
}
