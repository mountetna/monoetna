import { postRetrieveFiles, postCreateFolder } from '../api/files_api';

const addFiles = (files) => ({ type: 'ADD_FILES', files });
const addFolders = (folders) => ({ type: 'ADD_FOLDERS', folders });

export const retrieveFiles = ({folder_name}) => (dispatch) => {
  postRetrieveFiles(CONFIG.project_name, folder_name == undefined ? '' : folder_name)
    .then(({files, folders})=>{
      dispatch(addFiles(files));
      dispatch(addFolders(folders));
    })
}

export const createFolder = ({folder_name, parent_folder}) => (dispatch) => {
  postCreateFolder(CONFIG.project_name, parent_folder ? `${parent_folder}/${folder_name}` : folder_name)
    .then(({folders}) => dispatch(addFolders(folders)))
    .catch(error=>alert('Could not create folder!'));
}
