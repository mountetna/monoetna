import {postAuthorizeUpload} from '../api/upload_api';
import {errorMessage, message as showDialog} from './message_actions';
import {
  AddUploadCommand,
  CancelUploadCommand,
  UnpauseUploadCommand,
  Upload,
  PauseUploadCommand
} from '../workers/uploader';

export const ADD_FILES = 'ADD_FILES';
export const WORK = 'WORK';
export const WORK_FAILED = 'WORK_FAILED';
export const UNPAUSE_UPLOAD = 'UNPAUSE_UPLOAD';
export {PAUSE_UPLOAD, CANCEL_UPLOAD} from '../workers/uploader';

export const dispatchUploadWork = (dispatch, command) =>
  dispatch({
    type: WORK,
    work_type: 'upload',
    command
  });

export const fileSelected = ({file, folder_name, bucket_name}) => (
  dispatch,
  getState
) => {
  // Note: webkitRelativePath is NOT on the standards track, although it is supported in Firefox and Chrome currently.
  // There doesn't currently (June 2020) seem to be a standards way of getting these relative paths of directory uploads
  // other than using this attribute, which seems to work well enough for our target audiences for now.
  // TODO: Check for more standards way of doing this in the future, or if IE support is required.
  const relativePath = file.webkitRelativePath;
  const dest = relativePath || file.name;
  const {project_name} = CONFIG;
  let file_name = [folder_name, dest].filter((_) => _).join('/');

  // In the case of a directory upload, we need to first ensure that a containing directory exists.
  if (relativePath) {
    const parentDir = relativePath.split('/').slice(0, -1).join('/');
    if (parentDir) {
      return performMkdir(parentDir).then(performUpload);
    }
  }

  return performUpload();

  function performMkdir(parentDir) {
    return dispatch({
      type: 'CREATE_FOLDER',
      bucket_name,
      folder_name: parentDir,
      parent_folder: folder_name
    });
  }

  function performUpload() {
    return postAuthorizeUpload(
      window.location.origin,
      project_name,
      bucket_name,
      file_name
    )
      .then(({url}) => {
        dispatchUploadWork(
          dispatch,
          AddUploadCommand(
            Upload({
              file_name,
              url,
              file,
              project_name
            })
          )
        );
      })
      .catch(
        errorMessage(dispatch, 'warning', 'Upload failed', (error) => error)
      )
      .catch(
        errorMessage(
          dispatch,
          'error',
          'Upload failed',
          (error) => `Something bad happened: ${error}`
        )
      );
  }
};

export const unpauseUpload = ({upload}) => (dispatch) => {
  dispatchUploadWork(dispatch, UnpauseUploadCommand(upload));
};

export const pauseUpload = ({upload}) => (dispatch) => {
  dispatchUploadWork(dispatch, PauseUploadCommand(upload));
};

export const uploadComplete = ({upload}) => (dispatch) => {
  let {file} = upload;
  dispatch({type: ADD_FILES, files: [file]});
};

export const cancelUpload = ({upload}) => (dispatch) => {
  if (
    upload.status !== 'complete' &&
    !confirm('Are you sure you want to remove this upload?')
  )
    return;
  dispatchUploadWork(dispatch, CancelUploadCommand(upload));
};

export const uploadError = ({message_type, message, title, upload}) => (
  dispatch
) => {
  dispatch(showDialog(message_type, title, message));
};
