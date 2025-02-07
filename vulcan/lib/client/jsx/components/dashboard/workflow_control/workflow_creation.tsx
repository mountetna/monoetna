import React, {useCallback, useState, useContext} from 'react';

import Tooltip from '@material-ui/core/Tooltip';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';
import DialogActions from '@material-ui/core/DialogActions';
import AddIcon from '@material-ui/icons/Add';
import SaveIcon from '@material-ui/icons/Save';
import CancelIcon from '@material-ui/icons/Cancel';
import {VulcanContext} from '../../../contexts/vulcan_context';
import TextField from '@material-ui/core/TextField';
import Grid from '@material-ui/core/Grid';
import { updateWorkflowsWorkspaces } from '../../../actions/vulcan_actions';

const useStyles = makeStyles((theme) => ({
  dialog: {
    maxWidth: '90vw',
    maxHeight: '90vh'
  },
  subtitle: {display: 'inline'},
  button: {
    margin: '1rem'
  },
  loading: {
    marginLeft: '1rem'
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

export default function WorkflowCreateButtonModal({projectName}: {
  projectName: string;
}) {
  const classes = useStyles();
  const [repoUrl, setRepoUrl] = useState('');
  const [workflowName, setWorkflowName] = useState('');
  const [workflowNameError, setWorkflowNameError] = useState(false);
  const [open, setOpen] = useState(false);

  let {state,
    dispatch,
    showErrors,
    getWorkflows,
    createWorkflow
  } = useContext(VulcanContext);

  const handleCreate = useCallback((workflowName: string, repoUrl: string) => {
    // Check if workflow name is available
    const current_workflows = showErrors(getWorkflows(projectName));
    const current_names = Object.values(current_workflows).map(val => val.name)
    if (current_names.includes(workflowName)) {
      setWorkflowNameError(true);
    } else {
      showErrors(
        createWorkflow(
          projectName, repoUrl, workflowName
        )
      );
      setOpen(false);
      setWorkflowName('');
      setRepoUrl('');
      setWorkflowNameError(false);
      dispatch(updateWorkflowsWorkspaces());
    };
  }, [projectName, getWorkflows, createWorkflow]);

  function handleClose() {
    setOpen(false);
  }

  return (
    <>
      <Tooltip title='Add Workflow'>
        <Button
          className={classes.button}
          onClick={() => setOpen(true)}
          startIcon={<AddIcon/>}
          color='primary'
          variant='contained'
        >
          Add Workflow
        </Button>
      </Tooltip>
      <Dialog
          open={open}
          onClose={handleClose}
          maxWidth='xl'
        >
        <DialogTitle>
          Add a Workflow to this project
        </DialogTitle>
        <DialogContent className={classes.dialog}>
          <Grid container direction='column'>
            <Grid item>
              <Typography className={classes.helpdoc}>Please designate a valid github repo, and give the workflow a name.</Typography>
            </Grid>
            <Grid item>
              <TextField
                value={workflowName}
                label="Workflow Name"
                error={workflowNameError}
                helperText={workflowNameError ? 'Workflow name already exists for this project' : ''}
                InputLabelProps={{ shrink: true }}
                onChange={(event) => setWorkflowName(event.target.value)}
                size="small"
              />
            </Grid>
            <Grid item>
              <TextField
                value={repoUrl}
                label="Workflow URL"
                placeholder="github.com/<org>/<repo>"
                InputLabelProps={{ shrink: true }}
                onChange={(event) => setRepoUrl(event.target.value)}
                size="small"
              />
            </Grid>
          </Grid>
        </DialogContent>
        <DialogActions>
          <Tooltip title='Create Workflow'>
            <Button
              className={classes.propagateButton}
              onClick={() => handleCreate(workflowName, repoUrl)}
              startIcon={<SaveIcon/>}
              color='primary'
              variant='contained'
            >
              Create Workflow
            </Button>
          </Tooltip>
          <Tooltip title='Cancel'>
            <Button
              className={classes.propagateButton}
              onClick={handleClose}
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

