import { uploadAuthorize } from '../api/upload_api';

/*
 * Call to get approval to make an action on Metis.
 */
export const authorizeFile = ({ uploadFile }) => (dispatch) => {
  let { projectName, fileName, key } = uploadFile;

  if (!projectName || !fileName) {
    alert('The data to upload is not complete.');
    return;
  }

  uploadAuthorize(uploadFile)
    .then( response => response.text())
    .then( url => {
      // first set the upload url
      dispatch({ type: 'FILE_UPLOAD_AUTHORIZED', url, key });

      // then tell the worker to initialize the file
      dispatch({ type: 'WORK', worker: 'init', uploadFile, url });
    })
    .catch(
      () => alert('The upload could not be authorized.')
    )
}
