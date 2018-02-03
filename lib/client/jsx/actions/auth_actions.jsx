import { uploadAuthorize } from '../api/upload_api';

/*
 * Call to get approval to make an action on Metis.
 */
// the actual action
  
export const authorizeFile = (action) => (dispatch) => {
  let { uploadFile } = action;

  let requestItems = [ 'projectName', 'fileName' ];

  let uploadValid = requestItems.every(item => (uploadFile[item] != undefined));

  if (!uploadValid) {
    alert('The data to upload is not complete.');
    return;
  }

  uploadAuthorize(uploadFile)
    .then(
      response => dispatch(
        { type: 'FILE_UPLOAD_AUTHORIZED', response }
      )
    )
    .catch(
      () => alert('The upload could not be authorized.')
    )
}

export const fileUploadAuthorized = (action) => (dispatch) => {
  this.initializeFile(action.response.request);
}
