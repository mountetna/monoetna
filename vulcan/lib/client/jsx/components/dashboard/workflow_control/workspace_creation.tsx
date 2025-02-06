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
import {CreateWorkspaceResponse, Workflow, Workspace} from '../../../api_types';
import { VulcanContext } from '../../../contexts/vulcan_context';
import { runPromise, useAsyncCallback } from 'etna-js/utils/cancellable_helpers';
import { VulcanState } from '../../../reducers/vulcan_reducer';
import { includesClassNamePredicate } from '../../../test_utils/rendered';
import { workspaceId, workflowId } from '../../../selectors/workflow_selectors';
import Autocomplete from '@material-ui/lab/Autocomplete';

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
  workflow,
  workspaces,
}: {
  projectName: string;
  workflow: VulcanState['workflow'];
  workspaces: VulcanState['workspaces'];
}) {
  if (!workflow || !workflow.id) return null;
  const classes = useStyles();
  const invoke = useActionInvoker();
  let {state,
    showErrors,
    createWorkspace
  } = useContext(VulcanContext);

  const [workspaceName, setWorkspaceName] = useState('');
  const [branch, setBranch] = useState('main');
  const [repoVersion, setRepoVersion] = useState('')
  const [open, setOpen] = useState(false);

  const [handleCreateWorkspace] = useAsyncCallback(function* () {
    if (!workflow || !workflow.id) return;
    // ToDo: UI needed for branch / commit choice!
    const newSession: CreateWorkspaceResponse = yield* runPromise(showErrors(createWorkspace(projectName, workflow.id, workspaceName, branch, repoVersion)))
    invoke(
      pushLocation(
        `/${projectName}/workspace/${newSession.workspace_id}`
      )
    );
  }, [invoke, projectName, workflow, workspaceName, branch, repoVersion]);

  const defaultName = useMemo(() => {
    return workflow ? workflow.name : ''
  }, [workflow])

  // const pastVersions = useMemo(() => {
  //   return [
  //     ...new Set(workspaces
  //     .filter((w: Workspace) => {w.workflow_id == workflow.id})
  //     .map((w: Workspace) => w.git_version))
  //   ];
  // }, [workflow, workspaces])

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
          <Grid container direction='column' spacing={2}>
            <Grid item>
              <Typography className={classes.helpdoc}>This will establish a working directory on the remote compute server, and lauch a browser session where you will be able to run the workflow.</Typography>
            </Grid>
            <Grid item>
              <TextField
                value={workspaceName}
                label='Workspace Name'
                helperText="You will be able to change this later."
                InputLabelProps={{ shrink: true }}
                onChange={(event) => setWorkspaceName(event.target.value)}
                size="small"
              />
            </Grid>
            <Grid item>
              {/* ToDo: This could be a drop-down based on options given by the back-end */}
              <TextField
                value={branch}
                multiline
                label='Workflow branch'
                helperText='Ex: master or main'
                error={branch===''}
                InputLabelProps={{ shrink: true }}
                onChange={(event) => setBranch(event.target.value)}
                size="small"
              />
            </Grid>
            <Grid item>
              {/* ToDo: This could be a drop-down based on options given by the back-end */}
              <TextField
                label='Workflow Version'
                helperText='Commit SHA or Tag'
                error={repoVersion===''}
                InputLabelProps={{ shrink: true }}
                onChange={(event) => setRepoVersion(event.target.value)}
                size="small"
              />
              {/* <Autocomplete
                fullWidth
                freeSolo
                value={repoVersion}
                options={pastVersions}
                renderInput={(params: any) => (
                  <TextField
                    label='Workflow Version'
                    helperText='Commit SHA or Tag'
                    error={repoVersion===''}
                    InputLabelProps={{shrink: true}}
                    variant='outlined'
                    size="small"
                  />
                )}
                filterOptions={(options: string[], state: any) => {
                  let regex = new RegExp(state.inputValue);
                  return options.filter((o) => regex.test(o));
                }}
                onChange={(e: any, v: string | null) => {setRepoVersion(v==null ? '' : v)}}
              /> */}
            </Grid>
          </Grid>
        </DialogContent>
        <DialogActions>
          <Tooltip title='Create Workspace'>
            <Button
              className={classes.propagateButton}
              onClick={handleCreateWorkspace}
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
