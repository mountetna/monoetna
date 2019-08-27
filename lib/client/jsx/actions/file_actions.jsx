import {
  postRetrieveFiles, postProtectFile, postUnprotectFile, postRenameFile, deleteFile
} from '../api/files_api';

const addFiles = (files) => ({ type: 'ADD_FILES', files });
const addFolders = (folders) => ({ type: 'ADD_FOLDERS', folders });
const removeFiles = (files) => ({ type: 'REMOVE_FILES', files });

export const retrieveFiles = ({folder_name, bucket_name}) => (dispatch) =>
  postRetrieveFiles(CONFIG.project_name, bucket_name || 'files', folder_name == undefined ? '' : folder_name)
    .then(({files, folders})=>{
      dispatch(addFiles(files));
      dispatch(addFolders(folders));
    })
    .catch(error => dispatch({ type: 'INVALID_FOLDER' }));

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

