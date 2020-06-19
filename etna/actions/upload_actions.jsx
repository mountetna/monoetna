import { selectUpload, selectUploads } from '../selectors/directory-selector';
import { postAuthorizeUpload } from '../api/upload_api';
import { errorMessage } from '../actions/message_actions';

export const ADD_UPLOAD = 'ADD_UPLOAD';
export const ADD_FILES = 'ADD_FILES';
export const UPLOAD_STATUS = 'UPLOAD_STATUS';
export const REMOVE_UPLOAD = 'REMOVE_UPLOAD';
export const WORK = 'WORK';
export const UNQUEUE_UPLOADS = 'UNQUEUE_UPLOADS';
export const CANCEL_UPLOAD = 'CANCEL_UPLOAD';
export const CONTINUE_UPLOAD = 'CONTINUE_UPLOAD';
export const PAUSE_UPLOAD = 'PAUSE_UPLOAD';
export const WORK_FAILED = 'WORK_FAILED';


const work = (dispatch, command, args) => dispatch({
  type: WORK, work_type: 'upload', command,
  ...args
});

const getUpload = (state, file_name) => selectUpload(
  state, { project_name: CONFIG.project_name, file_name }
);

export const fileSelected = ({ file, folder_name, bucket_name }) => (dispatch, getState) => {
  // Note: webkitRelativePath is NOT on the standards track, although it is supported in Firefox and Chrome currently.
  // There doesn't currently (June 2020) seem to be a standards way of getting these relative paths of directory uploads
  // other than using this attribute, which seems to work well enough for our target audiences for now.
  // TODO: Check for more standards way of doing this in the future, or if IE support is required.
  const relativePath = file.webkitRelativePath;
  const dest = relativePath || file.name;
  let file_name = [ folder_name, dest ].filter(_=>_).join('/');

  // In the case of a directory upload, we need to first ensure that a containing directory exists.
  if (relativePath) {
    const parentDir = relativePath.split("/").slice(0, -1).join('/');
    if (parentDir) {
      return performMkdir(parentDir).then(performUpload);
    }
  }

  return performUpload();

  function performMkdir(parentDir) {
    return dispatch({ type: 'CREATE_FOLDER', bucket_name, folder_name: parentDir, parent_folder: folder_name });
  }

  function performUpload() {
    return postAuthorizeUpload(window.location.origin, CONFIG.project_name, bucket_name, file_name)
      .then(
        ({url}) => {
          dispatch({type: ADD_UPLOAD, project_name: CONFIG.project_name, file, file_name, url});
          let upload = getUpload(getState(), file_name);

          work(dispatch, 'start', {upload});
        }
      )
      .catch(
        errorMessage(dispatch, 'warning', 'Upload failed', error => error)
      )
      .catch(
        errorMessage(dispatch, 'error', 'Upload failed',
          error => `Something bad happened: ${error}`)
      );
  }
}

export const uploadStarted = ({ file_name }) => (dispatch, getState) => {
  let upload = getUpload(getState(), file_name);

  if (upload.status == 'active') work(dispatch, 'continue', { upload });
}

export const uploadBlobCompleted = ({file_name}) => (dispatch, getState) => {
  let upload = getUpload(getState(), file_name);

  if (upload.status == 'active') work(dispatch, 'continue', { upload });
}

export const uploadFileCompleted = ({upload}) => (dispatch) => {
  let { file } = upload;
  dispatch({ type: ADD_FILES, files: [ file ] });

  dispatch({ type: UNQUEUE_UPLOADS });
}

export const continueUpload = ({upload}) => (dispatch) => {
  dispatch({ type: UPLOAD_STATUS, upload, status: 'active' });
  work(dispatch, 'continue', { upload });
};

export const pauseUpload = ({upload}) => (dispatch) => {
  dispatch({ type: UPLOAD_STATUS, upload, status: 'paused' });
}

export const cancelUpload = ({upload}) => (dispatch) => {
  if (upload.status == 'complete') {
    dispatch({ type: REMOVE_UPLOAD, upload });
    return;
  }

  if (!confirm('Are you sure you want to remove this upload?')) return;

  dispatch({ type: WORK, work_type: 'upload', command: 'cancel', upload });
}

export const uploadFileCanceled = ({upload}) => (dispatch) => {
  dispatch({ type: REMOVE_UPLOAD, upload });
  dispatch({ type: UNQUEUE_UPLOADS });
}

export const unqueueUploads = () => (dispatch, getState) => {
  let uploads = Object.values(selectUploads(getState()));
  work(dispatch, 'unqueue', { uploads });
}
