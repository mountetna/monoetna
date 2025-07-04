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

import {VulcanContext} from '../../contexts/vulcan_context';
import {
  clearRunTriggers,
  setAttemptingToRun,
  setRunning,
  setWorkspace,
} from '../../actions/vulcan_actions';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
// import Vignette from '../vignette';
import VulcanHelp from './drawers/vulcan_help';
import { workflowName, workspaceFromRaw } from '../../selectors/workflow_selectors';
import {useWorkspace} from '../../contexts/workspace_context';
// import {json_get} from 'etna-js/utils/fetch';
import useUserHooks from '../../contexts/useUserHooks';
import Tag from '../dashboard/tag';
import Grid from '@material-ui/core/Grid';
import { useDataSync, useRunSyncing } from './data_sync';

// import RevisionHistory from 'etna-js/components/revision-history';

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

export default function WorkspaceManager() {
  const [saving, setSaving] = useState(false);
  const {
    state,
    dispatch,
    showError,
    showErrors,
    updateWorkspace,
    postUIValues,
    getFileNames,
    readFiles,
    getWorkspace,
    pullRunStatus,
    getIsRunning,
    requestRun,
    // updateFigure,
    // createFigure,
    // clearLocalSession
  } = useContext(VulcanContext);
  
  const {workflow, workspace, workspaceId, hasPendingEdits, complete} = useWorkspace();
  const {canEdit, guest} = useUserHooks();
  console.log({state});

  // const [modalIsOpen, setIsOpen] = useState(false);
  const [vulcanHelpIsOpen, setVulcanHelpIsOpen] = useState(false);
  const {workQueueable: committedStepPending, projectName, configId, isRunning} = state;

  const [localTags, setLocalTags] = useState<string[]>(workspace.tags || []);
  const [openTagEditor, setOpenTagEditor] = useState(false);
  // const [openRevisions, setOpenRevisions] = useState<boolean | null>(null);
  const [localTitle, setLocalTitle] = useState(workspace.name);
  const [updating, setUpdating] = useState(false);
  const [updatingTags, setUpdatingTags] = useState(false);
  const [updatingTitle, setUpdatingTitle] = useState(false);

  useEffect(() => {
    setLocalTitle(workspace.name);
  }, [workspace.name])

  function resetTags() {
    setLocalTags(workspace.tags || []);
  }
  useEffect(() => {
    resetTags();
  }, [workspace.tags])

  useDataSync(state, dispatch, showError, showErrors, getFileNames, readFiles, postUIValues);
  const {
    requestRunPolling,
    cancelRunning
  } = useRunSyncing(projectName, workspaceId, state.runId, showError, showErrors, pullRunStatus, getIsRunning, dispatch);

  const classes = useStyles();

  const workflow_name = workflowName(workflow);
  // const openModal = useCallback(() => setIsOpen(true), [setIsOpen]);
  // const closeModal = useCallback(() => setIsOpen(false), [setIsOpen]);;

  const run = useCallback(() => {
    if (!workspaceId || !configId) {
      showError('Possible UI Bug?: Missing info needed for requesting workspace run')
      return
    };
    if (isRunning) {
      showError('Possible UI Bug?: Cannot request run initiation while already running')
      return
    };
    dispatch(setAttemptingToRun(true));
    showErrors(requestRun(projectName, workspaceId, configId), (e) => {dispatch(setAttemptingToRun(false))})
    .then(runResponse => {
      dispatch(setRunning(runResponse.run_id));
    })
  }, [workspaceId, configId])
  useEffect(() => {
    if (state.isRunning) {
      showErrors(requestRunPolling());
    }
  }, [state.isRunning, state.isSyncing, requestRunPolling])

  const stop = useCallback(() => {
    // ToDo: Make and hook up a cancel_running api!
    // This does nothing currently
    cancelRunning();
  }, [cancelRunning]);

  // ToDo Later: once we figure out revisions.
  // const cancelSaving = useCallback(() => {
  //   setSaving(false);
  // }, []);

  // // ToDo Later: Save is only viable once we are handling revisions. (& Tag / title editing handled separately)
  // const handleSave = useCallback(
  //   (figure: VulcanFigure, newTags: string[], newTitle: string | undefined) => {
  //     let params = {
  //       workspaceId ,
  //       workflow_name: workflow_name,
  //       inputs: {...session.inputs},
  //       title: newTitle,
  //       tags: [...newTags]
  //     };

  //     if (!params.title) {
  //       params.title = prompt('Set a title for this figure') || undefined;
  //       if (!params.title) return;
  //     }

  //     if (params.figure_id) {
  //       params.comment = prompt('Enter a revision comment') || undefined;
  //       if (!params.comment) return;
  //     }

  //     setSaving(true);
  //     if (params.figure_id) {
  //       showErrors(
  //         updateFigure(session.project_name, params)
  //           .then((figureResponse) => {
  //             dispatch(setSessionAndFigure(figureResponse));
  //           })
  //           .finally(cancelSaving)
  //       );
  //     }
  //   },
  //   [
  //     workflow_name,
  //     session,
  //     cancelSaving,
  //     showErrors,
  //     updateFigure,
  //     invoke,
  //     dispatch,
  //     clearLocalSession,
  //     createFigure
  //   ]
  // );

  // const saveSession = useCallback(() => {
  //   if (hasPendingEdits) {
  //     if (!confirm('Pending edits will be discarded when saving. Proceed?')) {
  //       return;
  //     }
  //   }

  //   handleSave(figure, tags, localTitle);
  // }, [hasPendingEdits, handleSave, figure, tags, localTitle]);

  // // ToDo Later: Awaiting copyWorkspace API!
  // const copyFigure = useCallback(() => {
  //   let clone = {
  //     ...figure
  //   };

  //   delete clone.figure_id;

  //   handleSave(clone, [], `${localTitle} - copy`);
  // }, [figure, handleSave, localTitle]);

  // // ToDo Later: once we are handling revisions.
  // const loadRevision = useCallback(
  //   ({inputs, title, tags, id, workflow_snapshot}: VulcanRevision) => {
  //     dispatch(
  //       setSession({
  //         ...session,
  //         inputs,
  //         reference_figure_id: id
  //       } as VulcanFigureSession)
  //     );
  //     setLocalTitle(title);
  //     setTags(tags || []);
  //     requestPoll();
  //     setOpenRevisions(false);
  //     if (workflow_snapshot)
  //       {dispatch(setWorkflow(workflow_snapshot, session.project_name));}
  //   },
  //   [dispatch, session, requestPoll]
  // );

  const setAllNotUpdating = () => {
    setUpdating(false);
    setUpdatingTags(false);
    setUpdatingTitle(false);
  }

  const handleUpdateWorkspace = useCallback((title?: string, tags?: string[]) => {
    // ToDo: Handle this update from the api return instead, to be sure we stay accurate
    setUpdating(true);
    showErrors(updateWorkspace(projectName, workspaceId as number, title, tags), (e) => setAllNotUpdating())
      .then(() => {
        showErrors(getWorkspace(projectName, workspaceId as number), (e) => setAllNotUpdating())
        .then((workspaceRaw) => {
          dispatch(setWorkspace(workspaceFromRaw(workspaceRaw), projectName))
          setAllNotUpdating()
        });
      });
  }, [projectName, workspaceId])

  const handleCloseEditTags = useCallback(() => {
    resetTags();
    setOpenTagEditor(false);
  }, []);

  const isPublic = useMemo(() => (workspace.tags || []).includes('published'), [workspace.tags]);
  function togglePublicTag() {
    const currentTags = workspace.tags || [];
    const newTags = currentTags.includes('published') ?
    currentTags.filter((v)=>v!='published') :
      [...currentTags, 'published']
    setLocalTags(newTags);
    setUpdatingTags(true);
    handleUpdateWorkspace(undefined, newTags);
  }

  const running = useMemo(() => state.isRunning, [state.isRunning]);
  const disableRunButton =
    complete || running || (hasPendingEdits && !committedStepPending) || state.isSyncing || state.attemptingToRun;
  const disableRunReason = running ?
    'Workspace is already runnning' :
    state.attemptingToRun ?
    'Requesting to Run' :
    hasPendingEdits && !committedStepPending ?
    'An input is has pending edits' :
    state.isSyncing ?
    'Awaiting sync with remote workspace' :
    'no remaining work to run'

  // Catch auto-pass 'Run' trigger
  useEffect(() => {
    if (state.triggerRun.length > 0 && !isRunning) {
      dispatch(clearRunTriggers(state.triggerRun));
      run();
    }
  }, [state.triggerRun, dispatch, run]);

  // // ToDo Later: once we figure out revisions.
  // const inputsChanged = useMemo(() => {
  //   return !_.isEqual(paramValuesToRaw(status.params), status.last_params) ||
  //     !_.isEqual(status.ui_contents, pick(uiContentsFromFiles(workspace, status.file_contents), Object.keys(status.ui_contents)));
  // }, [status, workspace]);

  // // ToDo Later: once we allow updateWorkspace.
  // const titleChanged = useMemo(() => {
  //   return localTitle !== workspace.name;
  // }, [workspace.name, localTitle]);

  // // ToDo Later: once we allow updateWorkspace.
  // const tagsChanged = useMemo(() => {
  //   return !_.isEqual(tags, workspace.tags);
  // }, [workspace.tags, tags]);

  // // ToDo Later: once we are handling revisions.
  // const viewingRevision = useMemo(() => {
  //   return figure.id !== session.reference_figure_id;
  // }, [figure, session]);

  const canSave = false
  // const canSave = useMemo(() => {
  //   return (
  //     (titleChanged || inputsChanged || tagsChanged) &&
  //     !(running || saving) &&
  //     !viewingRevision
  //   );
  // }, [
  //   running,
  //   saving,
  //   inputsChanged,
  //   titleChanged,
  //   tagsChanged,
  //   viewingRevision
  // ]);

  const editor = useMemo(() => canEdit(workspace), [
    workspace,
    canEdit
  ]);

  if (!workflow_name || !workspace) return null;

  return (
    <div className='session-manager'>
      <div className='session-header'>
        <Breadcrumbs
          className='session-workflow-name'
          classes={{
            li: classes.title
          }}
        >
          <Link href={`/${workflow.project_name}`}>{workflow.project_name}</Link>
          <Typography>{workflow.name}</Typography>
          {editor ? (
            <Grid container direction='row'>
              <Grid item>
                <TextField
                  fullWidth
                  value={localTitle}
                  margin='none'
                  InputProps={{
                    disableUnderline: true,
                    inputProps: {
                      className: classes.titleText
                    }
                  }}
                  disabled={updating}
                  variant='standard'
                  onChange={(e) => setLocalTitle(e.target.value)}
                  placeholder='Untitled'
                />
              </Grid>
              {localTitle != workspace.name && <Grid item xs={2}>
                <FlatButton
                  className='header-btn-name-save'
                  icon={updatingTitle? 'spinner fa-spin' : 'save'}
                  label='Save'
                  title='Save Workspace Title'
                  disabled={updating}
                  onClick={() => {
                    setUpdatingTitle(true);
                    handleUpdateWorkspace(localTitle, undefined)
                  }}
                />
              </Grid>}
            </Grid>
          ) : (
            <Typography>{workspace.name}</Typography>
          )
          }
        </Breadcrumbs>
        <React.Fragment>
          <FlatButton
            icon='book'
            className='header-btn vignette'
            label='Vulcan'
            onClick={() => setVulcanHelpIsOpen(true)}
          />
          <ReactModal
            isOpen={vulcanHelpIsOpen}
            onRequestClose={() => setVulcanHelpIsOpen(false)}
            style={modalStyles}
            contentLabel='Vignette'
          >
            <VulcanHelp/>
          </ReactModal>
        </React.Fragment>
        {/* {workflow.vignette && (
          <React.Fragment>
            <FlatButton
              icon='book'
              className='header-btn vignette'
              label='Vignette'
              onClick={openModal}
            />
            <ReactModal
              isOpen={modalIsOpen}
              onRequestClose={closeModal}
              style={modalStyles}
              contentLabel='Vignette'
            >
              <Vignette workflowName={workflow_name} />
            </ReactModal>
          </React.Fragment>
        )} */}
        {running ? (
          <FlatButton
            className={'header-btn'}
            icon='stop'
            label='Stop'
            title='Coming soon, Cancel running work'
            onClick={stop}
            disabled={true}
          />
        ) : (
          <FlatButton
            className={'header-btn run'}
            icon={state.attemptingToRun ? 'spinner fa-spin' : 'play'}
            label='Run'
            title={disableRunButton ? disableRunReason : 'Run workflow'}
            onClick={() => run()}
            disabled={disableRunButton}
          />
        )}
        {editor ?
          <>
            <FlatButton
              className='header-btn save'
              icon='save'
              label='Save'
              title='Coming by Beta, create a Save Point you can return to later'
              // onClick={saveSession}
              disabled={!canSave || updating}
            />
            <FlatButton
              className='header-btn public-private'
              icon={updatingTags ? 'spinner fa-spin' : isPublic ? 'fa-solid fa-eye-slash' : 'fa-solid fa-eye'}
              label={`${isPublic ? 'Unpublish' : 'Publish'}`}
              title={`Make the current figure ${
                isPublic ? 'private to you' : 'public to all with project access'
              }`}
              onClick={() => {
                togglePublicTag();
              }}
              disabled={updating || guest(projectName || '')}
            />
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
                  disabled={_.isEqual(localTags,workspace.tags)}
                >
                  Save Tags
                </Button>
                <Button onClick={handleCloseEditTags} color='primary'>
                  Close
                </Button>
              </DialogActions>
            </Dialog>
          </> : <FlatButton
            className='header-btn copy'
            icon='copy'
            label='Copy'
            title='Copy current workflow parameters to new figure'
            disabled={true}
            onClick={() => {}}
          />
        }
        {/* editor ? (
          <>
            <FlatButton
              className='header-btn edit-tags'
              icon='history'
              label='Revisions'
              title='Revisions'
              onClick={() => setOpenRevisions(true)}
            />
            {openRevisions != null && (
              <RevisionHistory
                open={openRevisions}
                onClose={() => setOpenRevisions(false)}
                revisionDoc={({
                  inputs,
                  title,
                  tags,
                  dependencies
                }: VulcanRevision) =>
                  JSON.stringify({inputs, title, tags, dependencies}, null, 2)
                }
                update={loadRevision}
                getRevisions={() =>
                  json_get(
                    `/api/${workspace.project}/figure/${figure.figure_id}/revisions`
                  )
                }
              />
            )}
            // ToDo Later: If choose to allow workspaces to pull from a different commit, or updated container
            <AdvancedSessionControls session={session} figure={figure} />
          </>
        ) : (
          <FlatButton
            className='header-btn copy'
            icon='copy'
            label='Copy'
            title='Copy current workflow parameters to new figure'
            onClick={copyFigure}
          />
        )*/}
      </div>
      <div className='session-feed-container'>
        <InputFeed />
        <OutputFeed />
      </div>
    </div>
  );
}
