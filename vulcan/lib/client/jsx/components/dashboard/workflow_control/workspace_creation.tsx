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
import Switch from '@material-ui/core/Switch';
import InputLabel from '@material-ui/core/InputLabel';

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
  const [valUse, setValUse] = useState({version: '...awaiting...', lastUsed: 'never'});
  const [requestBy, setRequestBy] = useState<'branch' | 'tagOrSha'>('branch');
  const [open, setOpen] = useState(false);
  const [creating, setCreating] = useState(false);
  const [versionText, setVersionText] = useState('main');
  const [versionHelperText, setVersionHelperText] = useState('Branch name (e.g. main)');
  const [createTag, setCreateTag] = useState('Create Workspace');

  const [handleCreateWorkspace] = useAsyncCallback(function* (name: string, version: string, request_by: 'branch' | 'tagOsSha') {
    if (!workflow || !workflow.id) return;
    setCreating(true);
    showErrors(createWorkspace(projectName, workflow.id, name, version, request_by))
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
  }, [invoke, projectName, workflow, workspaceName]);

  const defaultName = useMemo(() => {
    return workflow ? workflow.name : ''
  }, [workflow])

  //ToDo: grab past tags & split out tag as requestBy option?
  const pastShas: {version: string, lastUsed: string}[] = useMemo(() => {
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
  const pastBranches: {version: string, lastUsed: string}[] = useMemo(() => {
    const wspaces: WorkspaceMinimal[] = workspaces.filter((w: WorkspaceMinimal) => w.workflow_id == workflow.id);
    const branches = [...new Set(wspaces
      .filter((w: WorkspaceMinimal) => w.git_request_by=='branch')
      .map((w: WorkspaceMinimal) => w.git_request) as string[]
    )];
    return branches.map((version) => {
      const used = wspaces.filter((w: WorkspaceMinimal) => w.git_request_by=='branch' && w.git_request == version)
        .map((w: WorkspaceMinimal) => w.created_at.split(' ')[0])
      return {
        version: version,
        // ToDo: better selection of latest!
        lastUsed: used[used.length-1]
      }
    })
  }, [workflow, workspaces]);
  useEffect(() => {
    // Find 'main' if exists for last used element
    if (requestBy == 'branch' && valUse.version == '...awaiting...' && pastBranches.length !== undefined) {
      const pastMain = pastBranches.filter(b => b.version=='main')
      const newBranch = {version: 'main', lastUsed: pastMain.length>0 ? pastMain[0].lastUsed : 'never'}
      setValUse(newBranch);
      setCreateTag('Create Workspace');
    }
  }, [requestBy, valUse, pastBranches])

  const helperBase = requestBy == 'branch' ? 'Branch name (e.g. main)' : 'Tag or Commit SHA';
  const pastUse = requestBy == 'branch' ? pastBranches : pastShas;
  const optsUse = [...pastUse, {version: requestBy == 'branch' ? '...awaiting...' : '', lastUsed: 'never'}]
  const versionUI = <Autocomplete
    key={`git-version-request by-${requestBy}`}
    freeSolo
    value={valUse}
    options={optsUse}
    inputValue={versionText}
    onInputChange={(event: any, value: string) => {
      setVersionText(value)
    }}
    renderInput={(params: any) => (
      <TextField
        {...params}
        label='Workflow Version'
        helperText={versionHelperText}
        error={valUse.version == '...awaiting...' || valUse.version==='' || valUse.version != versionText}
        InputLabelProps={{shrink: true}}
        size="small"
      />
    )}
    getOptionLabel={(option) => option.version}
    renderOption={
      (option: typeof valUse, state: object) => {
      return <Tooltip title={'last used: ' + option.lastUsed} placement='right'>
        <Typography>
          {option.version}
        </Typography>
      </Tooltip>
    }}
    filterOptions={(options: typeof pastShas, state: any) => {
      let regex = new RegExp(state.inputValue);
      return options.filter((o) => regex.test(o.version) && o.version!='' && o.version!='...awaiting...');
    }}
    onChange={(e: any, v: {version: string, lastUsed: string} | string | null) => {
      if (v != null && typeof v === 'object') {
        setVersionHelperText(helperBase);
        setValUse(v);
      };
      if (typeof v === 'string' && v != '') {
        const match = pastUse.filter(p => p.version==v);
        if (match.length>0) {
          setValUse(match[0]);
        } else {
          setValUse({version: v, lastUsed: 'never'});
        }
      }
    }}
  />

  useEffect(() => {
    if (valUse.version === '' || valUse.version == '...awaiting...') {
      setCreateTag('Workflow Version not set')
      setVersionHelperText(helperBase + ', *Hit Enter/Return to use the current value*')
    } else if (valUse.version != versionText) {
      setCreateTag('Workflow Version input in error sate')
      setVersionHelperText(helperBase + ', *Hit Enter/Return to use the current value*')
    } else if (createTag!='Create Workspace') {
      setVersionHelperText(helperBase);
      setCreateTag('Create Workspace');
    }
  }, [valUse, versionText, createTag])
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
              <InputLabel htmlFor='request-by-switch' shrink>Request Workspace Version By</InputLabel>
              <Typography component="div" key='request-by-switch'>
                <Grid component="label" container alignItems="center" spacing={1}>
                  <Grid item>Branch</Grid>
                  <Grid item>
                    <Switch
                      checked={requestBy=='tagOrSha'}
                      onChange={() => {
                        setValUse({version: requestBy != 'branch' ? '...awaiting...' : '', lastUsed: 'never'});
                        setVersionText(requestBy != 'branch' ? 'main' : '');
                        setRequestBy(requestBy!='branch' ? 'branch' : 'tagOrSha');
                      }}
                    />
                  </Grid>
                  <Grid item>TagOrSha</Grid>
                </Grid>
              </Typography>
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
                if (!disableCreate) handleCreateWorkspace(workspaceName, valUse.version, requestBy)
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
