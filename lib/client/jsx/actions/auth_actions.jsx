import { postAuthorizeUpload } from '../api/upload_api';

/*
 * Call to get approval to make an action on Metis.
 */

export const fileSelected = ({ file }) => (dispatch) => {
  postAuthorizeUpload(CONFIG.project_name, file.name)
    .then( response => response.text())
    .then( url => {
      // first set the upload url
      dispatch({ type: 'FILE_UPLOAD_AUTHORIZED', file, url });

      // then tell the worker to initialize the file
      dispatch({ type: 'WORK', worker: 'upload', command: 'init', file, url });
    })
    .catch(
      () => alert('The upload could not be authorized.')
    )
}
