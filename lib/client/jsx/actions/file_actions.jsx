import { postRetrieveFiles } from '../api/files_api';

export const retrieveFiles = () => (dispatch) => {
  postRetrieveFiles(CONFIG.project_name)
    .then( ({files}) => {
      // first set the upload url
      dispatch({ type: 'ADD_FILES', files });
    })
}
