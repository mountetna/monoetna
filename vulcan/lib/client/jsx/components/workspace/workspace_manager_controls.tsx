import React, {
  useCallback,
  useEffect,
  useState,
  useContext,
  useMemo
} from 'react';
import * as _ from 'lodash';
import ReactModal from 'react-modal';
import FlatButton from 'etna-js/components/flat-button';

import Grid from '@material-ui/core/Grid';
import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import Autocomplete from '@material-ui/lab/Autocomplete';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
// import {pushLocation} from 'etna-js/actions/location_actions';

import Breadcrumbs from '@material-ui/core/Breadcrumbs';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Link from '@material-ui/core/Link';
import Tooltip from '@material-ui/core/Tooltip';

import {WorkspaceContext} from '../../contexts/workspace_context';
import VulcanHelp from './drawers/vulcan_help';
import { hasScheduledSteps, workflowName, workspaceFromRaw } from '../../selectors/workflow_selectors';
import useUserHooks from '../../contexts/useUserHooks';
import Tag from '../dashboard/tag';
import { useDataSync, useRunSyncing } from './data_sync';
import Vignette from './drawers/vignette';

const modalStyles = {
  content: {
    top: '50%',
    left: '50%',
    right: 'auto',
    bottom: 'auto',
    marginRight: '-50%',
    transform: 'translate(-50%, -50%)'
  }
};

const useStyles = makeStyles((theme) => ({
  title: {
    '&:last-child': {flex: '1 1 auto'}
  },
  titleText: {
    textOverflow: 'ellipsis',
    overflow: 'hidden'
  },
  tags: {
    padding: '12.5px !important'
  },
  editTags: {
    width: '600px'
  }
}));

const ModalButton = ({icon, title, label, modalLabel, disabled, children}) => {
  const [ showModal, setShowModal ] = React.useState(false);

  return <>
    <FlatButton
      icon={icon}
      className='header-btn'
      label={label}
      title={title}
      disabled={disabled}
      onClick={() => {setShowModal(true)}}
    />
    <ReactModal
      isOpen={showModal}
      onRequestClose={() => setShowModal(false)}
      style={modalStyles}
      contentLabel={modalLabel}
    >
    { children }
    </ReactModal>
  </>
}

const VignetteButton = () => {
  const { state: { file_contents } } = React.useContext(WorkspaceContext);

  const hasVignette = 'vignette.md' in file_contents;

  return <ModalButton
    icon='book'
    label='Workflow'
    title={hasVignette ? 'Workflow Readme' : 'Workflow Readme Unavailable'}
    disabled={!hasVignette}>
    <Vignette/>
  </ModalButton>
}

const HelpButton = () => {
  return <ModalButton
    icon='book'
    label='Vulcan'
    title='Vulcan Interface Overview'>
    doop
    <VulcanHelp/>
  </ModalButton>
}

const RunButton = () => {
  const { state: { isRunning } } = React.useContext(WorkspaceContext);
  const disableRunButton = false;
  return <>
    {isRunning ? (
      <FlatButton
        className={'header-btn'}
        icon='stop'
        label='Stop'
        title='Cancel running work'
        onClick={stop}
        disabled={true}
      />
    ) : (
      <FlatButton
        className={'header-btn run'}
        icon='play'
        label='Run'
        title={'Run workflow'}
        onClick={() => run()}
        disabled={disableRunButton}
      />
    )}
  </>
}

const TagButton = () => {
  const updatingTags = false;
  const updating = false;
  const { state: { workspace } } = React.useContext(WorkspaceContext);
  const [ openTagEditor, setOpenTagEditor ] = React.useState(false);
  const [localTags, setLocalTags] = React.useState<string[]>(workspace.tags || []);

  const classes = useStyles();

  const handleCloseEditTags = () => setOpenTagEditor(false);

  return <>
    <FlatButton
      className='header-btn edit-tags'
      icon={updatingTags ? 'spinner fa-spin' : 'tags'}
      label='Edit tags'
      title='Edit tags'
      disabled={updating}
      onClick={() => setOpenTagEditor(true)}
    />
    <Dialog
      maxWidth='md'
      open={openTagEditor}
      onClose={handleCloseEditTags}
    >
      <DialogTitle id='tag-editor'>Edit Tags</DialogTitle>
      <DialogContent className={classes.editTags}>
        <Autocomplete
          fullWidth
          multiple
          freeSolo
          className='figure-edit-tag-autocomplete'
          classes={{
            input: classes.tags
          }}
          defaultValue={localTags}
          id='figure-edit-tags-filter'
          options={localTags.filter((t) => t!='published' && t!='highlighted').concat('highlighted')}
          renderInput={(params: any) => (
            <TextField {...params} label='Tags' variant='outlined' />
          )}
          renderTags={(tags: string[], getTagProps: any) =>
            tags.map((tag, index) => (
              <Tag {...getTagProps({index})} label={tag} />
            ))
          }
          renderOption={(option: string, state: any) => (
            <span>{option}</span>
          )}
          filterOptions={(options: string[], state: any) => {
            let regex = new RegExp(state.inputValue);
            return options.filter((o) => regex.test(o));
          }}
          onChange={(e: any, v: string[]) => setLocalTags(v)}
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={() => {
          setUpdatingTags(true);
          handleUpdateWorkspace(undefined, localTags)
          handleCloseEditTags()
          }}
          color='primary'
          disabled={_.isEqual(localTags,workspace?.tags)}
        >
          Save Tags
        </Button>
        <Button onClick={handleCloseEditTags} color='primary'>
          Close
        </Button>
      </DialogActions>
    </Dialog>
  </>
}

const WorkspaceManagerControls = () => {
  const classes = useStyles();

  const { canEdit } = useUserHooks();
  const { state: { workspace } } = React.useContext(WorkspaceContext);

  const editor = canEdit(workspace);

  const updatingTags = false;
  const isPublic = false;

  return <>
    <VignetteButton/>
    <HelpButton/>
    <RunButton/>
    {editor && 
      <>
        <FlatButton
          className='header-btn public-private'
          icon={updatingTags ? 'spinner fa-spin' : isPublic ? 'fa-solid fa-eye-slash' : 'fa-solid fa-eye'}
          label={`${isPublic ? 'Unpublish' : 'Publish'}`}
          title={`Make the current figure ${
            isPublic ? 'private to you' : 'public to all with project access'
          }`}
          onClick={() => {
            //togglePublicTag();
          }}
        />
        <TagButton/>
      </>
    }
  </>
}

export default WorkspaceManagerControls;
