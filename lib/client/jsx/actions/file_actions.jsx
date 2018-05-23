import { postRetrieveFiles, postCreateFolder } from '../api/files_api';

const addFiles = (dispatch) => ({files}) => dispatch({ type: 'ADD_FILES', files });

export const retrieveFiles = ({folder_name}) => (dispatch) => {
  postRetrieveFiles(CONFIG.project_name, folder_name == undefined ? '' : folder_name)
    .then(addFiles(dispatch))
}

export const createFolder = ({folder_name, parent_folder}) => (dispatch) => {
  postCreateFolder(CONFIG.project_name, `${parent_folder}/${folder_name}`)
    .then(addFiles(dispatch))
    .catch(error=>alert('Could not create folder!'));
}
