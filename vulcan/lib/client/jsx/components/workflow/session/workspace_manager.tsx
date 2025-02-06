import React, {
  useCallback,
  useEffect,
  useState,
  useContext,
  useMemo
} from 'react';
import * as _ from 'lodash';
// import ReactModal from 'react-modal';
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

import {VulcanContext} from '../../../contexts/vulcan_context';
import {
  clearCommittedStepPending,
  clearRunTriggers,
  setWorkspace,
  removeSync,
  updateFiles,
} from '../../../actions/vulcan_actions';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
// import Vignette from '../vignette';
import { allFilesToBuffer, filesReturnToMultiFileContent, uiContentsFromFiles, workflowName } from '../../../selectors/workflow_selectors';
import {useWorkspace} from '../../../contexts/workspace_context';
// import {json_get} from 'etna-js/utils/fetch';
import useUserHooks from '../../useUserHooks';
import Tag from '../../dashboard/tag';
import Grid from '@material-ui/core/Grid';
import {runPromise, useAsync} from 'etna-js/utils/cancellable_helpers';
import { MultiFileContentResponse, WorkspaceStatus } from '../../../api_types';
import { VulcanState } from '../../../reducers/vulcan_reducer';

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
    showErrors,
    requestPoll,
    cancelPolling,
    updateWorkspace,
    getFileNames,
    readFiles,
    // updateFigure,
    // createFigure,
    // clearLocalSession
  } = useContext(VulcanContext);
  const {workflow, workspace, workspaceId, hasPendingEdits, complete} = useWorkspace();
  const {canEdit, guest} = useUserHooks();
  console.log({state});

  // const [modalIsOpen, setIsOpen] = useState(false);
  const {status, workQueueable: committedStepPending, projectName} = state;

  const [tags, setTags] = useState<string[]>(workspace.tags || []);
  const [openTagEditor, setOpenTagEditor] = useState(false);
  // const [openRevisions, setOpenRevisions] = useState<boolean | null>(null);
  const [localTitle, setLocalTitle] = useState(workspace.name);

  useEffect(() => {
    setLocalTitle(workspace.name);
  }, [workspace.name])

  const invoke = useActionInvoker();

  const classes = useStyles();

  const workflow_name = workflowName(workflow);
  // const openModal = useCallback(() => setIsOpen(true), [setIsOpen]);
  // const closeModal = useCallback(() => setIsOpen(false), [setIsOpen]);

  useAsync(function* () {
    if (status.to_sync.length > 0) {
      // Push files / params to compute server for  step
      const pushStep = [...status.to_sync][0]
      requestPoll(state, true, false, pushStep);
      dispatch(removeSync(pushStep));
    }
  }, [status.to_sync])

  useAsync(function* () {
    if (state.update_files && !!workspaceId) {
      // Push files / params to compute server for  step
      const fileNamesRaw: {files: string[]} = yield* runPromise(showErrors(getFileNames(projectName, workspaceId)));
      const fileNames = fileNamesRaw.files;
      const update = {
        output_files: fileNames,
        file_contents: {},
        ui_contents: {}
      }
      const fileGrabs = allFilesToBuffer(workspace).filter(f => fileNames.includes(f));
      if (fileGrabs.length > 0) {
        const filesContentRaw: MultiFileContentResponse = yield* runPromise(readFiles(projectName, workspaceId, fileGrabs))
        const filesContent = filesReturnToMultiFileContent(filesContentRaw);
        update['file_contents'] = filesContent;
        update['ui_contents'] = uiContentsFromFiles(workspace, filesContent);
      } else {
        update['file_contents'] = {};
        update['ui_contents'] = uiContentsFromFiles(workspace);
      }
      dispatch(updateFiles(update))
    }
  }, [state.update_files, workspaceId])

  const run = useCallback((_state: VulcanState) => {
    showErrors(requestPoll(_state,false,true));
    dispatch(clearCommittedStepPending());
    // ToDo: Additional accounting per completed steps!
  }, [requestPoll, dispatch, showErrors]);
  const stop = useCallback(() => {
    cancelPolling()
    // ToDo: Additional accounting per completed steps!
  }, [cancelPolling]);

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

  const handleUpdateWorkspace = useCallback((title?: string, tags?: string[]) => {
    // ToDo: Handle this update from the api return instead, to be sure we stay accurate
    const newWorkspace = {...workspace}
    if (!!title && !_.isEqual(title, workspace.name)) {
      newWorkspace['name'] = title;
    }
    if (!!tags && !_.isEqual(tags, workspace.tags)) {
      newWorkspace['tags'] = tags;
    }
    showErrors(
      updateWorkspace(projectName, workspaceId as number, title, tags)
        // .then((workspaceResponse) => {
        //   dispatch(setWorkspace(workspaceResponse, workspaceResponse.project));
        // })
        .then(() => {
          dispatch(setWorkspace(newWorkspace, projectName));
        })
    );
  }, [projectName, workspaceId])

  const handleCloseEditTags = useCallback(() => {
    setOpenTagEditor(false);
  }, []);

  function togglePublicTag() {
    if (tags.includes('public')) {
      setTags(tags.filter((v)=>v!='public'));
    } else {
      setTags([...tags, 'public']);
    }
  }

  const running = useMemo(() => state.pollingState > 0, [state.pollingState]);
  const disableRunButton =
    complete || running || (hasPendingEdits && !committedStepPending);

  // Catch auto-pass 'Run' trigger
  useEffect(() => {
    if (state.triggerRun.length > 0) {
      dispatch(clearRunTriggers(state.triggerRun));
      run(state);
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

  const isPublic = useMemo(() => (tags || []).includes('public'), [tags]);

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
                  variant='standard'
                  onChange={(e) => setLocalTitle(e.target.value)}
                  placeholder='Untitled'
                />
              </Grid>
              {localTitle != workspace.name && <Grid item xs={2}>
                <FlatButton
                  className='header-btn-name-save'
                  icon='save'
                  label='Save'
                  title='Save Workspace Title'
                  onClick={() => handleUpdateWorkspace(localTitle, undefined)}
                />
              </Grid>}
            </Grid>
          ) : (
            <Typography>{workspace.name}</Typography>
          )
          }
        </Breadcrumbs>
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
        {state.pollingState ? (
          <FlatButton
            className={'header-btn'}
            icon='stop'
            label='Stop'
            title='Stop workflow'
            onClick={stop}
          />
        ) : (
          <FlatButton
            className={'header-btn run'}
            icon='play'
            label='Run'
            title='Run workflow'
            onClick={() => run(state)}
            disabled={disableRunButton}
          />
        )}
        {editor ?
          <>
            <FlatButton
              className='header-btn public-private'
              icon={`${isPublic ? 'lock' : 'unlock'}`}
              label={`Make ${isPublic ? 'private' : 'public'}`}
              title={`Make the current figure ${
                isPublic ? 'private' : 'public'
              }`}
              onClick={() => {
                togglePublicTag();
                handleUpdateWorkspace(undefined, tags);
              }}
              disabled={guest(projectName || '')}
            />
            <FlatButton
              className='header-btn edit-tags'
              icon='tags'
              label='Edit tags'
              title='Edit tags'
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
                  defaultValue={tags}
                  id='figure-edit-tags-filter'
                  options={tags}
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
                  onChange={(e: any, v: string[]) => setTags(v)}
                />
              </DialogContent>
              <DialogActions>
                <Button onClick={() => handleUpdateWorkspace(undefined, tags)} color='primary'>
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
              className='header-btn save'
              icon='save'
              label='Save'
              title='Save current workflow parameters to current figure'
              onClick={saveSession}
              disabled={!canSave}
            />
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
