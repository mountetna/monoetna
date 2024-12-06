import React, {useCallback, useContext, useMemo, useState} from 'react';

import Tooltip from '@material-ui/core/Tooltip';
import Button from '@material-ui/core/Button';
import WorkIcon from '@material-ui/icons/Work';
import SaveIcon from '@material-ui/icons/Save';
import CancelIcon from '@material-ui/icons/Cancel';
import Typography from '@material-ui/core/Typography';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';
import DialogActions from '@material-ui/core/DialogActions';
import TextField from '@material-ui/core/TextField';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';

import {pushLocation} from 'etna-js/actions/location_actions';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {Workflow} from '../../../api_types';
import { VulcanContext } from '../../../contexts/vulcan_context';

const useStyles = makeStyles((theme) => ({
  dialog: {
    maxWidth: '90vw',
    maxHeight: '90vh'
  },
  subtitle: {display: 'inline'},
  button: {
    margin: '1rem'
  },
  helpdoc: {
    maxWidth: '600px',
    marginTop: '1rem',
    marginBottom: '1rem'
  },
  propagateButton: {
    marginBottom: '1rem'
  }
}));

export default function WorkspaceCreateButtonModal({
  projectName,
  workflow
}: {
  projectName: string;
  workflow: Workflow | null;
}) {
  if (!workflow || !workflow.id) return;
  const classes = useStyles();
  const invoke = useActionInvoker();
  let {state,
    showErrors,
    createWorkspace
  } = useContext(VulcanContext);

  const [branch, setBranch] = useState('main');
  const [workspaceName, setWorkspaceName] = useState('');
  const [open, setOpen] = useState(false);

  const handleCreateWorkspace = useCallback(() => {
    if (!workflow || !workflow.id) return;
    // ToDo: UI needed for branch / commit choice!
    const newSession = showErrors(createWorkspace(projectName, workflow.id, branch, workspaceName))
    invoke(
      pushLocation(
        `/${projectName}/${workflow.name}/${newSession.workspace_id}`
      )
    );
  }, [invoke, projectName, workflow]);

  const defaultName = useMemo(() => {
    return workflow ? workflow.name : ''
  }, [workflow])

  return (
    <>
      <Tooltip title='Create Workspace'>
        <Button
          disabled={!workflow}
          onClick={() => {
            setWorkspaceName(defaultName);
            setOpen(true);
          }}
          color='primary'
          startIcon={<WorkIcon />}
        >
          Create Workspace
        </Button>
      </Tooltip>
      <Dialog
          open={open}
          onClose={() => setOpen(false)}
          maxWidth='xl'
        >
        <DialogTitle>
          Create New {workflow.name} Workspace 
        </DialogTitle>
        <DialogContent className={classes.dialog}>
          <Grid container>
            <Grid item>
              <Typography className={classes.helpdoc}>This will establish a working directory on the remote compute server, and lauch a browser session where you will be able to run the workflow.</Typography>
            </Grid>
            <Grid item>
              <TextField
                value={workspaceName}
                multiline
                label="Workspace Name"
                helperText="You will be able to change this later."
                InputLabelProps={{ shrink: true }}
                onChange={(event) => setWorkspaceName(event.target.value)}
                size="small"
              />
            </Grid>
            <Grid item>
              {/* ToDo: This should be a drop-down based on options given by the back-end */}
              <TextField
                value={branch}
                multiline
                label="Workflow branch"
                error={branch===''}
                InputLabelProps={{ shrink: true }}
                onChange={(event) => setBranch(event.target.value)}
                size="small"
              />
            </Grid>
          </Grid>
        </DialogContent>
        <DialogActions>
          <Tooltip title='Create Workspace'>
            <Button
              className={classes.propagateButton}
              onClick={() => handleCreateWorkspace}
              startIcon={<SaveIcon/>}
              color='primary'
              variant='contained'
            >
              Create Workspace
            </Button>
          </Tooltip>
          <Tooltip title='Cancel'>
            <Button
              onClick={() => setOpen(false)}
              startIcon={<CancelIcon/>}
              color='secondary'
              variant='contained'
            >
              Cancel
            </Button>
          </Tooltip>
        </DialogActions>
      </Dialog>
    </>
  );
}
