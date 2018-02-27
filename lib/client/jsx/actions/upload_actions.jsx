import { postUploadBlob, postUploadStart } from '../api/upload_api';

export const queueUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'WORK', worker: 'upload', command: 'start', upload });
}

export const pauseUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'WORK', worker: 'upload', command: 'pause', upload });
}

export const cancelUpload = ({upload}) => (dispatch) => {
  if(!confirm('Are you sure you want to remove this upload?')) return;
  dispatch({ type: 'WORK', worker: 'upload', command: 'cancel', upload });
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
