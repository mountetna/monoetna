import React, {useCallback} from 'react';

import Tooltip from '@material-ui/core/Tooltip';
import Button from '@material-ui/core/Button';
import InsertChartIcon from '@material-ui/icons/InsertChart';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import {makeStyles} from '@material-ui/core/styles';

import {pushLocation} from 'etna-js/actions/location_actions';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {Workflow} from '../../api_types';

const useStyles = makeStyles((theme) => ({
  controls: {
  }
}));

export default function WorkflowControls({
  project_name,
  workflow
}: {
  project_name: string;
  workflow: Workflow | null;
}) {
  const classes = useStyles();
  const invoke = useActionInvoker();

  const handleCreateFigure = useCallback(() => {
    if (!workflow) return;
    invoke(
      pushLocation(
        `/${project_name}/figure/new/${workflow.name.replace('.cwl', '')}`
      )
    );
  }, [invoke, project_name, workflow]);

  return (
    <ButtonGroup
      variant='contained'
      color='primary'
      className={classes.controls}
    >
      <Tooltip title='Create figure'>
        <Button
          disabled={!workflow}
          onClick={handleCreateFigure}
          color='primary'
          startIcon={<InsertChartIcon />}
        >
          Create figure
        </Button>
      </Tooltip>
    </ButtonGroup>
  );
}
