import { selectUpload } from '../selectors/directory-selector';

const work = (dispatch, worker_name, command, args) => dispatch({
  type: 'WORK', work_type: 'upload', worker_name, command,
  ...args
});

const getUpload = (getState, file_name) => selectUpload(
  getState(), { project_name: CONFIG.project_name, file_name }
);

export const fileSelected = ({ file, folder_name, bucket_name }) => (dispatch) => {
  let file_name = [ folder_name, file.name ].filter(_=>_).join('/');

  work(dispatch, file_name, 'authorize', {
    base_url: window.location.origin,
    project_name: CONFIG.project_name, bucket_name, file, file_name
  });
}

export const uploadAuthorized = ({ file, file_name, url }) => (dispatch, getState) => {
  dispatch({ type: 'ADD_UPLOAD', file, file_name, url });

  let upload = getUpload(getState, file_name);

  work(dispatch, file_name, 'start', { upload });
}

export const uploadStarted = ({ file_name }) => (dispatch, getState) => {
  let upload = getUpload(getState, file_name);

  if (upload.status == 'active')
    work(dispatch, file_name, 'continue', { upload });
}

export const uploadBlobCompleted = ({file_name}) => (dispatch, getState) => {
  let upload = getUpload(getState, file_name);

  if (upload.status == 'active')
    work(dispatch, file_name, 'continue', { upload });
}

export const uploadFileCompleted = ({upload}) => (dispatch) => {
  let { file } = upload;
  dispatch({ type: 'ADD_FILES', files: [ file ] });
}

export const continueUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'UPLOAD_STATUS', upload, status: 'active' });
  work(dispatch, upload.file_name, 'continue', { upload });
}

export const pauseUpload = ({upload}) => (dispatch) => {
  dispatch({ type: 'UPLOAD_STATUS', upload, status: 'paused' });
}

export const cancelUpload = ({upload}) => (dispatch) => {
  if(upload.status != 'complete' && !confirm('Are you sure you want to remove this upload?')) return;
  dispatch({ type: 'WORK', work_type: 'upload', worker_name: upload.file_name, command: 'cancel', upload });
}
