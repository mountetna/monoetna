import React, {useContext, useEffect, useMemo, useState} from 'react';

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
import Autocomplete from '@material-ui/lab/Autocomplete';

import {pushLocation} from 'etna-js/actions/location_actions';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {WorkspaceMinimal} from '../../../api_types';
import { VulcanContext } from '../../../contexts/vulcan_context';
import { useAsyncCallback } from 'etna-js/utils/cancellable_helpers';
import { VulcanState } from '../../../reducers/vulcan_reducer';
import LoadingIcon from '../loading_icon';

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
  workflow: VulcanState['workflow'] | null;
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
  const [repoVersion, setRepoVersion] = useState({version: '', lastUsed: 'never'})
  const [open, setOpen] = useState(false);
  const [creating, setCreating] = useState(false);
  const [versionText, setVersionText] = useState(repoVersion.version)
  const [versionHelperText, setVersionHelperText] = useState('Commit SHA')
  const [createTag, setCreateTag] = useState('Create Workspace')

  const [handleCreateWorkspace] = useAsyncCallback(function* () {
    if (!workflow || !workflow.id) return;
    setCreating(true);
    showErrors(createWorkspace(projectName, workflow.id, workspaceName, branch, repoVersion.version))
    .then((newSession) => {
      setCreating(false);
      invoke(
        pushLocation(
          `/${projectName}/workspace/${newSession.workspace_id}`
        )
      );
    })
    .catch((e) => {
      setCreating(false)
      // ToDo: Make use of error to highlight when the user gets their branch or version wrong.
      // if (e[0].includes(`branch ${branch} not found`)) {
      //   setBranchError(true)
      // } else if (e[0].includes('did not match any file(s) known to git')) {
      //   setVersionError(true)
      // }
    })
  }, [invoke, projectName, workflow, workspaceName, branch, repoVersion]);

  const defaultName = useMemo(() => {
    return workflow ? workflow.name : ''
  }, [workflow])

  const pastVersions: {version: string, lastUsed: string}[] = useMemo(() => {
    const wspaces: WorkspaceMinimal[] = workspaces.filter((w: WorkspaceMinimal) => w.workflow_id == workflow.id);
    const versions = [...new Set(wspaces.map((w: WorkspaceMinimal) => w.git_version) as string[])]
    return versions.map((version) => {
      const used = wspaces.filter((w: WorkspaceMinimal) => w.git_version == version)
        .map((w: WorkspaceMinimal) => w.created_at.split(' ')[0])
      return {
        version: version,
        // ToDo: better selection of latest!
        lastUsed: used[used.length-1]
      }
    })
  }, [workflow, workspaces])

  const versionUI = <Autocomplete
    freeSolo
    value={repoVersion}
    options={[...pastVersions, {version: '', lastUsed: 'never'}]}
    inputValue={versionText}
    onInputChange={(event: any, value: string) => {
      if (value != repoVersion.version) {
        setVersionHelperText('Commit SHA, *Hit Enter/Return to use the current value*')
      } else {
        setVersionHelperText('Commit SHA')
      }
      setVersionText(value)
    }}
    renderInput={(params: any) => (
      <TextField
        {...params}
        label='Workflow Version'
        helperText={versionHelperText}
        error={repoVersion.version==='' || repoVersion.version != versionText}
        InputLabelProps={{shrink: true}}
        size="small"
      />
    )}
    getOptionLabel={(option) => option.version}
    renderOption={
      (option: typeof repoVersion, state: object) => {
      return <Tooltip title={'last used: ' + option.lastUsed} placement='right'>
        <Typography>
          {option.version}
        </Typography>
      </Tooltip>
    }}
    filterOptions={(options: typeof pastVersions, state: any) => {
      let regex = new RegExp(state.inputValue);
      return options.filter((o) => regex.test(o.version) && o.version!='');
    }}
    onChange={(e: any, v: {version: string, lastUsed: string} | string | null) => {
      if (v != null && typeof v === 'object') {
        setVersionHelperText('Commit SHA');
        setRepoVersion(v);
      };
      if (typeof v === 'string' && v != '') {
        const match = pastVersions.filter(p => p.version==v);
        if (match.length>0) {
          setRepoVersion(match[0]);
        } else {
          setRepoVersion({version: v, lastUsed: 'never'});
        }
      }
    }}
  />

  useEffect(() => {
    if (repoVersion.version === '') {
      setCreateTag('Workflow Version not set')
    } else if (repoVersion.version != versionText) {
      setCreateTag('Workflow Version input in error sate')
    } else if (createTag!='Create Workspace') {
      setCreateTag('Create Workspace')
    }
  }, [repoVersion, versionText, createTag])
  const disableCreate = createTag!='Create Workspace';

  return (
    <>
      <Tooltip title='Create Workspace'>
        <Button
          className={classes.button}
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
                fullWidth
                label='Workspace Name'
                helperText="You will be able to change this later."
                InputLabelProps={{ shrink: true }}
                onChange={(event) => setWorkspaceName(event.target.value)}
                size="small"
              />
            </Grid>
            <Grid item>
              {versionUI}
            </Grid>
          </Grid>
        </DialogContent>
        <DialogActions>
          <Tooltip title={createTag} placement='top'>
            <Button
              className={classes.propagateButton}
              onClick={() => {
                if (!disableCreate) handleCreateWorkspace()
              }}
              startIcon={creating ? <LoadingIcon/> : <SaveIcon/>}
              color={disableCreate ? 'secondary' : 'primary'}
              variant='contained'
            >
              Create Workspace
            </Button>
          </Tooltip>
          <Tooltip title='Cancel' placement='top'>
            <Button
              className={classes.propagateButton}
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
