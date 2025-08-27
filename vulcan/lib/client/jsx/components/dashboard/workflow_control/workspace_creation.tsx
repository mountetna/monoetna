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
import FlatButton from 'etna-js/components/flat-button';

import git from "isomorphic-git";
import http from "isomorphic-git/http/web";

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
  const classes = useStyles();
  const invoke = useActionInvoker();
  const {
    showErrors,
    createWorkspace
  } = useContext(VulcanContext);

  const [workspaceName, setWorkspaceName] = useState('');
  const [valUse, setValUse] = useState({version: '...awaiting...', lastUsed: 'never'});
  const [branches, setBranches] = useState<string[]>([]);
  const [tags, setTags] = useState<string[]>([]);
  const [requestBy, setRequestBy] = useState<'branch' | 'tagOrSha'>('branch');
  const [open, setOpen] = useState(false);
  const [creating, setCreating] = useState(false);
  const [versionText, setVersionText] = useState('main');
  const [versionHelperText, setVersionHelperText] = useState('Branch name (e.g. main)');
  const [createTag, setCreateTag] = useState('Create Workspace');
  const [advanced, setAdvanced] = useState(false);

  const [handleCreateWorkspace] = useAsyncCallback(function* (name: string, version: string) {
    if (!workflow || !workflow.id) return;
    setCreating(true);
    showErrors(createWorkspace(projectName, workflow.id, name, version), (e) => {setCreating(false); setOpen(false)})
    .then((newSession) => {
      setCreating(false);
      invoke(
        pushLocation(
          `/${projectName}/workspace/${newSession.workspace_id}`
        )
      );
    });
  }, [invoke, projectName, workflow, workspaceName]);

  const defaultName = useMemo(() => {
    return workflow ? workflow.name : ''
  }, [workflow])

  useEffect(() => {
    if (!workflow) return;
    showErrors(
      git.getRemoteInfo({ http, url: workflow.repo_remote_url, corsProxy: 'https://cors.isomorphic-git.org' }),
      (e) => {
        setBranches(['main'])
        setTags([])
      }
    )
    .then((remote: any) => {
      setBranches(Object.keys(remote['refs']['heads']))
      setTags( ('tags' in remote['refs']) ? Object.keys(remote['refs']['tags']) : []);
    })
  }, [workflow])

  // Options to show
  type VersionOpts = {
    branches: {version: string, lastUsed: string}[];
    tags_or_commits: {version: string, lastUsed: string}[];
  };
  const versionOpts: VersionOpts = useMemo(() => {
    const out: VersionOpts = {branches: [], tags_or_commits: []}
    if (!workflow) return out;
    const wspaces: WorkspaceMinimal[] = workspaces
      .filter((w: WorkspaceMinimal) => w.workflow_id == workflow.id)
      .sort((a,b) => a.created_at < b.created_at ? 1 : -1);
    const versions: string[] = [...new Set(wspaces.map((w: WorkspaceMinimal) => w.git_ref) as string[])]
    const unusedBranches = [...branches];
    const unusedTags = [...tags];
    for (let ind in versions) {
      let v = versions[ind];
      let k: 'branches' | 'tags_or_commits' = 'tags_or_commits'
      if (branches.includes(v)) {
        k = 'branches';
        unusedBranches.splice(unusedBranches.indexOf(v), 1);
      }
      if (tags.includes(v)) {
        unusedTags.splice(unusedTags.indexOf(v), 1);
      }
      let used = wspaces
        .filter((w: WorkspaceMinimal) => w.git_ref == v)
        .map((w: WorkspaceMinimal) => w.created_at)
        // .sort((a,b) => a < b ? 1 : -1)
        [0].split(' +')[0]
      out[k].push({version: v, lastUsed: used})
    };
    if (unusedBranches.length>0) {
      for (let ind in unusedBranches) {
        out['branches'].push({version: unusedBranches[ind], lastUsed: 'never'})
      }
    }
    if (unusedTags.length>0) {
      for (let ind in unusedTags) {
        out['tags_or_commits'].push({version: unusedTags[ind], lastUsed: 'never'})
      }
    }
    return out;
  }, [workflow, workspaces, branches, tags])

  useEffect(() => {
    // Find 'main' if exists for last used element
    if (valUse.version == '...awaiting...' && requestBy == 'branch' && versionOpts['branches'].length > 0) {
      const pastMain = versionOpts['branches'].filter(b => b.version=='main')
      const newBranch = {version: 'main', lastUsed: pastMain.length>0 ? pastMain[0].lastUsed : 'never'}
      setValUse(newBranch);
      setCreateTag('Create Workspace');
    }
  }, [requestBy, valUse, versionOpts])

  // Creation button label & disabling, and update versionUI helper text
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

  if (!workflow || !workflow.id) {
    return <Tooltip title='Select a Workflow to enable'>
        <Button
          className={classes.button}
          onClick={() => {}}
          size='large'
          startIcon={<WorkIcon />}
          color='secondary'
        >
          Create New Workspace
        </Button>
      </Tooltip>
  };

  // versionUI
  const helperBase = requestBy == 'branch' ? 'Branch name (e.g. main)' : 'Tag or Commit SHA';
  const pastUse = requestBy == 'branch' ? versionOpts['branches'] : versionOpts['tags_or_commits'];
  const optsUse = [...pastUse, {version: requestBy == 'branch' ? '...awaiting...' : '', lastUsed: 'never'}]
  const optionDisplay = (option: typeof valUse, state: object) => {
    return <Grid container spacing={3}>
      <Grid item>
        <Typography>
          {option.version}
        </Typography>
      </Grid>
      <Grid item>
        <Typography color='secondary' variant='body2'>
          {`Last Requested: ${option.lastUsed}`}
        </Typography>
      </Grid>
    </Grid>
  };
  const onSelect = (e: any, v: typeof valUse | string | null) => {
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
  };
  const advancedVersionUI = <Grid container style={{paddingLeft: '20px'}}>
    <Grid item xs={5}>
      <InputLabel htmlFor='request-by-switch' disabled={!advanced} shrink>
        Git Reference Type
      </InputLabel>
      <Typography component="div" key='request-by-switch' variant='body2' color={advanced ? 'textPrimary' : 'textSecondary'}>
        <Grid component="label" container alignItems="center">
          <Grid item>Branch</Grid>
          <Grid item>
            <Switch
              checked={requestBy=='tagOrSha'}
              color="default"
              disabled={!advanced}
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
    <Grid item xs={7}>
      <Autocomplete
        key={`git-version-request by-${requestBy}`}
        freeSolo
        fullWidth
        value={valUse}
        disabled={!advanced}
        options={optsUse}
        inputValue={versionText}
        onInputChange={(event: any, value: string) => {
          setVersionText(value)
        }}
        renderInput={(params: any) => (
          <TextField
            {...params}
            label='Git Reference'
            helperText={versionHelperText}
            error={valUse.version == '...awaiting...' || valUse.version==='' || valUse.version != versionText}
            InputLabelProps={{shrink: true}}
            size="small"
          />
        )}
        getOptionLabel={(option) => {
          if (typeof option === 'string') return option
          return option.version
        }}
        renderOption={optionDisplay}
        filterOptions={(options: typeof optsUse, state: any) => {
          let regex = new RegExp(state.inputValue);
          return options.filter((o) => regex.test(o.version) && !['','...awaiting...'].includes(o.version))
        }}
        onChange={onSelect}
      />
    </Grid>
  </Grid>
  const disableCreate = createTag!='Create Workspace' || creating;

  return (
    <>
      <Tooltip title={`Create New \'${workflow.name}\' Workspace`}>
        <Button
          className={classes.button}
          disabled={!workflow}
          onClick={() => {
            setWorkspaceName(defaultName);
            setOpen(true);
          }}
          size='large'
          color='primary'
          startIcon={<WorkIcon />}
        >
          Create New Workspace
        </Button>
      </Tooltip>
      <Dialog
          open={open}
          onClose={() => setOpen(false)}
          maxWidth='xl'
        >
        <DialogTitle>
          {`Create New \'${workflow.name}\' Workspace`} 
        </DialogTitle>
        <DialogContent className={classes.dialog}>
          <Grid container direction='column' spacing={4}>
            <Grid item>
              <Typography className={classes.helpdoc} variant='body1'>
                This will set up a "workspace", a folder on our compute server where you can run the workflow.
                After it's ready, we'll send you over to your Workspace Control Panel.
              </Typography>
              <Typography className={classes.helpdoc} variant='body2' style={{marginTop: '1rem'}}>
                The choices below are optional.
              </Typography>
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
              <InputLabel htmlFor='version-basic-advanced' shrink>
                Workflow Version
              </InputLabel>
              <Grid
                key='version-basic-advanced'
                container alignItems='center'
                justifyContent='flex-start'
              >
                <Grid item>
                  <Typography variant='body1' color={advanced ? 'textSecondary' : 'textPrimary'}>
                    Default
                  </Typography>
                </Grid>
                <Grid item>
                  <Switch
                    checked={advanced}
                    onChange={() => {
                      setValUse({version: '...awaiting...' , lastUsed: 'never'});
                      setVersionText('main');
                      setRequestBy('branch');
                      setAdvanced(!advanced)
                    }}
                  />
                </Grid>
                <Grid item>
                  <Typography variant='body2' color={advanced ? 'textPrimary' : 'textSecondary'}>
                    Advanced
                  </Typography>
                </Grid>
              </Grid>
              {advancedVersionUI}
            </Grid>
          </Grid>
        </DialogContent>
        <DialogActions>
          <Tooltip title={createTag} placement='top'>
            <Button
              className={classes.propagateButton}
              onClick={() => {
                if (!disableCreate) handleCreateWorkspace(workspaceName, valUse.version)
              }}
              startIcon={creating ? <LoadingIcon/> : <SaveIcon/>}
              color={disableCreate ? 'secondary' : 'primary'}
              variant='contained'
            >
              <strong>Create Workspace</strong>
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
