import { postAuthorizeUpload } from '../api/upload_api';

/*
 * Call to get approval to make an action on Metis.
 */
export const authorizeUpload = ({ upload }) => (dispatch) => {
  let { project_name, file_name, file, key } = upload;

  if (!project_name || !file_name) {
    alert('The data to upload is not complete.');
    return;
  }

  postAuthorizeUpload(upload)
    .then( response => response.text())
    .then( url => {
      // first set the upload url
      dispatch({ type: 'FILE_UPLOAD_AUTHORIZED', url, key });

      // then tell the worker to initialize the file
      dispatch({ type: 'WORK', worker: 'init', file, url });
    })
    .catch(
      () => alert('The upload could not be authorized.')
    )
}
