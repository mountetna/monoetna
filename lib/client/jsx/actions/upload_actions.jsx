import { postUploadBlob, postUploadStart } from '../api/upload_api';

export const queueUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'WORK', worker: 'upload', command: 'start', upload });
}

export const pauseUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'WORK', worker: 'upload', command: 'pause', upload });
}

export const cancelUpload = () => (dispatch) => {
  if(!confirm('Are you sure you want to remove this upload?')) return;
  this.uploadWorker.postMessage({ command: 'cancel' });
}

export const recoverUpload = () => (dispatch) => {
  this.recoverUpload(action.uploadFile, action.fileMetadata);
}

export const startUpload = ({ upload, url }) => (dispatch) => {
  postUploadStart(url, upload)
    .then( response =>
      dispatch({ type: 'FILE_UPLOAD_STATUS', upload: response, status: 'paused' })
    )
    .catch(
      () => alert('The upload could not be started.')
    )
}
