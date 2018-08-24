import { postRetrieveFiles, postCreateFolder, postProtectFile, postUnprotectFile, postRenameFile, deleteFile } from '../api/files_api';

const addFiles = (files) => ({ type: 'ADD_FILES', files });
const removeFiles = (files) => ({ type: 'REMOVE_FILES', files });
const addFolders = (folders) => ({ type: 'ADD_FOLDERS', folders });

export const retrieveFiles = ({folder_name}) => (dispatch) =>
  postRetrieveFiles(CONFIG.project_name, folder_name == undefined ? '' : folder_name)
    .then(({files, folders})=>{
      dispatch(addFiles(files));
      dispatch(addFolders(folders));
    })
    .catch(error => dispatch({ type: 'INVALID_FOLDER' }));

export const createFolder = ({folder_name, parent_folder}) => (dispatch) =>
  postCreateFolder(CONFIG.project_name, parent_folder ? `${parent_folder}/${folder_name}` : folder_name)
    .then(({folders}) => dispatch(addFolders(folders)))
    .catch(error=>alert('Could not create folder!'));

export const removeFile = ({file}) => (dispatch) => {
  if (!confirm(`Are you sure you want to remove ${file.file_path}?`)) return;

  deleteFile(
    CONFIG.project_name, file.file_path
  )
    .then(({files}) => dispatch(removeFiles(files)))
    .catch(error=>alert('Could not remove file!'));
}

export const protectFile = ({file}) => (dispatch) => {
  postProtectFile(
    CONFIG.project_name, file.file_path
  )
    .then(({files}) => dispatch(addFiles(files)))
    .catch(error=>alert('Could not protect file!'));
}

export const unprotectFile = ({file}) => (dispatch) => {
  if (!confirm(`Are you sure you want to unprotect ${file.file_path}?`)) return;

  postUnprotectFile(
    CONFIG.project_name, file.file_path
  )
    .then(({files}) => dispatch(addFiles(files)))
    .catch(error=>alert('Could not unprotect file!'));
}

export const renameFile = ({file, new_file_path}) => (dispatch) => {
  postRenameFile(
    CONFIG.project_name, file.file_path, new_file_path
  )
    .then(({files}) => {
      dispatch(removeFiles([file]));
      dispatch(addFiles(files));
    })
    .catch(error=>alert('Could not rename file!'));
}
