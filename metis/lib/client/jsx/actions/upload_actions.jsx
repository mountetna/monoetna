import { selectUpload, selectUploads } from '../selectors/directory-selector';
import { postAuthorizeUpload } from '../api/upload_api';
import { errorMessage } from '../actions/message_actions';

const work = (dispatch, command, args) => dispatch({
  type: 'WORK', work_type: 'upload', command,
  ...args
});

const getUpload = (state, file_name) => selectUpload(
  state, { project_name: CONFIG.project_name, file_name }
);

export const fileSelected = ({ file, folder_name, bucket_name }) => (dispatch, getState) => {
  let file_name = [ folder_name, file.name ].filter(_=>_).join('/');

  postAuthorizeUpload(window.location.origin, CONFIG.project_name, bucket_name, file_name)
    .then(
      ({url}) => {
        dispatch({ type: 'ADD_UPLOAD', file, file_name, url });

        let upload = getUpload(getState(), file_name);

        work(dispatch, 'start', { upload });
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
  dispatch({ type: 'ADD_FILES', files: [ file ] });

  dispatch({ type: 'UNQUEUE_UPLOADS' });
}

export const continueUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'UPLOAD_STATUS', upload, status: 'active' });
  work(dispatch, 'continue', { upload });
}

export const pauseUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'UPLOAD_STATUS', upload, status: 'paused' });
}

export const cancelUpload = ({upload}) => (dispatch) => {
  if (upload.status == 'complete') {
    dispatch({ type: 'REMOVE_UPLOAD', upload });
    return;
  }

  if (!confirm('Are you sure you want to remove this upload?')) return;

  dispatch({ type: 'WORK', work_type: 'upload', command: 'cancel', upload });
}

export const uploadFileCanceled = ({upload}) => (dispatch) => {
  dispatch({ type: 'REMOVE_UPLOAD', upload});
  dispatch({ type: 'UNQUEUE_UPLOADS' });
}

export const unqueueUploads = () => (dispatch, getState) => {
  let uploads = Object.values(selectUploads(getState()));
  work(dispatch, 'unqueue', { uploads });
}
