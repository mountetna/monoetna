import {
  postCreateFolder, postProtectFolder, postUnprotectFolder, postRenameFolder, deleteFolder
} from '../api/folders_api';

const addFolders = (folders) => ({ type: 'ADD_FOLDERS', folders });
const removeFolders = (folders) => ({ type: 'REMOVE_FOLDERS', folders });

export const createFolder = ({bucket_name, folder_name, parent_folder}) => (dispatch) =>
  postCreateFolder(CONFIG.project_name, bucket_name, parent_folder ? `${parent_folder}/${folder_name}` : folder_name)
    .then(({folders}) => dispatch(addFolders(folders)))
    .catch(error=>alert('Could not create folder!'));

export const removeFolder = ({bucket_name, folder}) => (dispatch) => {
  if (!confirm(`Are you sure you want to remove ${folder.folder_path}?`)) return;

  deleteFolder(
    CONFIG.project_name, bucket_name, folder.folder_path
  )
    .then(({folders}) => dispatch(removeFolders(folders)))
    .catch(error=>alert('Could not remove folder!'));
}

export const protectFolder = ({bucket_name, folder}) => (dispatch) => {
  postProtectFolder(
    CONFIG.project_name, bucket_name, folder.folder_path
  )
    .then(({folders}) => dispatch(addFolders(folders)))
    .catch(error=>alert('Could not protect folder!'));
}

export const unprotectFolder = ({bucket_name, folder}) => (dispatch) => {
  if (!confirm(`Are you sure you want to unprotect ${folder.folder_path}?`)) return;

  postUnprotectFolder(
    CONFIG.project_name, bucket_name, folder.folder_path
  )
    .then(({folders}) => dispatch(addFolders(folders)))
    .catch(error=>alert('Could not unprotect folder!'));
}

export const renameFolder = ({bucket_name, folder, new_folder_path}) => (dispatch) => {
  postRenameFolder(
    CONFIG.project_name, bucket_name, folder.folder_path, new_folder_path
  )
    .then(({folders}) => {
      dispatch(removeFolders([folder]));
      dispatch(addFolders(folders));
    })
    .catch(error=>alert('Could not rename folder!'));
}
