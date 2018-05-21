import { postRetrieveFiles, postCreateFolder } from '../api/files_api';

export const retrieveFiles = ({folder_name}) => (dispatch) => {
  postRetrieveFiles(CONFIG.project_name, folder_name == undefined ? '' : folder_name)
    .then( ({files}) => {
      // first set the upload url
      dispatch({ type: 'ADD_FILES', files });
    })
}

export const createFolder = ({folder_name}) => (dispatch) => {
  postCreateFolder(CONFIG.project_name, folder_name)
    .then( ({folders}) => {
      dispatch({ type: 'ADD_FOLDERS', folders })
    })
}
