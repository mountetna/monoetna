import {
  postCreateFolder, postProtectFolder, postUnprotectFolder, postRenameFolder, deleteFolder,
  getTouchFolder
} from '../api/folders_api';
import { errorMessage } from './message_actions';

const addFolders = ({folders, bucket_name, current_folder}) => ({ type: 'ADD_FOLDERS', folders, bucket_name, current_folder });
const removeFolders = (folders) => ({ type: 'REMOVE_FOLDERS', folders });

export const createFolder = ({bucket_name, folder_name, parent_folder}) => (dispatch) =>
  postCreateFolder(CONFIG.project_name, bucket_name, parent_folder ? `${parent_folder}/${folder_name}` : folder_name)
    .then(({folders}) => dispatch(addFolders({folders})))
    .catch(
      errorMessage(dispatch, 'warning', 'Folder creation failed', error => error)
    );

export const removeFolder = ({bucket_name, folder, recursive}) => (dispatch) => {
  deleteFolder(
    CONFIG.project_name, bucket_name, folder.folder_path, recursive
  )
    .then(({folders}) => dispatch(removeFolders(folders)))
    .catch(
      response => response.then(
        ({error}) => {
          if (error == 'Folder is not empty') {
            if (!confirm(`"${folder.folder_path}" is not empty. Are you sure you want to remove the folder and its contents?`)) return;
            removeFolder({bucket_name, folder, recursive: true})(dispatch);
            return;
          }
          errorMessage(dispatch, 'warning', 'Folder removal failed', error => error)(error);
        }
      )
    );
}


export const protectFolder = ({bucket_name, folder}) => (dispatch) => {
  postProtectFolder(
    CONFIG.project_name, bucket_name, folder.folder_path
  )
    .then(({folders}) => dispatch(addFolders({folders})))
    .catch(
      errorMessage(dispatch, 'warning', 'Folder protection failed', error => error)
    );
}

export const unprotectFolder = ({bucket_name, folder}) => (dispatch) => {
  if (!confirm(`Are you sure you want to unprotect ${folder.folder_path}?`)) return;

  postUnprotectFolder(
    CONFIG.project_name, bucket_name, folder.folder_path
  )
    .then(({folders}) => dispatch(addFolders({folders})))
    .catch(
      errorMessage(dispatch, 'warning', 'Folder unprotection failed', error => error)
    );
}

export const renameFolder = ({bucket_name, folder, new_folder_path, new_bucket_name, current_folder}) => (dispatch) => {
  postRenameFolder(
    CONFIG.project_name, bucket_name, folder.folder_path, new_folder_path, new_bucket_name
  )
    .then(({folders}) => {
      dispatch(removeFolders([folder]));
      dispatch(addFolders({folders, bucket_name, current_folder}));
    })
    .catch(
      errorMessage(dispatch, 'warning', 'Folder renaming failed', error => error)
    );
}

export const touchFolder = ({bucket_name, folder}) => (dispatch) => {
  getTouchFolder(
    CONFIG.project_name, bucket_name, folder.folder_path
  )
    .then(({folders}) => {
      dispatch(removeFolders([folder]));
      dispatch(addFolders({folders}));
    })
    .catch(
      errorMessage(dispatch, 'warning', 'Folder touching failed', error => error)
    );
}

