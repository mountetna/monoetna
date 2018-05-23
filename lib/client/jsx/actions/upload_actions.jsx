import { postAuthorizeUpload } from '../api/upload_api';

/*
 * Call to get approval to make an action on Metis.
 */

export const fileSelected = ({ file, folder_name }) => (dispatch) => {
  let file_name = [ folder_name, file.name ].filter(_=>_).join('/');
  postAuthorizeUpload(CONFIG.project_name, file_name)
    .then( response => response.text())
    .then( url => {
      // first set the upload url
      dispatch({ type: 'FILE_UPLOAD_AUTHORIZED', file, file_name, url });

      // then tell the worker to initialize the file
      dispatch({ type: 'WORK', worker: 'upload', command: 'init', file, url });
    })
    .catch(
      () => alert('The upload could not be authorized.')
    )
}

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
