import {
  postRetrieveBuckets, postRetrieveFiles, postProtectFile, postUnprotectFile, postRenameFile, deleteFile,
  postCreateFolder, postProtectFolder, postUnprotectFolder, postRenameFolder, deleteFolder,
} from '../api/files_api';

const addFiles = (files) => ({ type: 'ADD_FILES', files });
const addBuckets = (buckets) => ({ type: 'ADD_BUCKETS', buckets });
const removeFiles = (files) => ({ type: 'REMOVE_FILES', files });
const addFolders = (folders) => ({ type: 'ADD_FOLDERS', folders });
const removeFolders = (folders) => ({ type: 'REMOVE_FOLDERS', folders });

export const retrieveFiles = ({folder_name, bucket_name}) => (dispatch) =>
  postRetrieveFiles(CONFIG.project_name, bucket_name || 'files', folder_name == undefined ? '' : folder_name)
    .then(({files, folders})=>{
      dispatch(addFiles(files));
      dispatch(addFolders(folders));
    })
    .catch(error => dispatch({ type: 'INVALID_FOLDER' }));

export const retrieveBuckets = () => (dispatch) =>
  postRetrieveBuckets(CONFIG.project_name).then(
    ({buckets}) => dispatch(addBuckets(buckets))
  );

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

export const createFolder = ({folder_name, parent_folder}) => (dispatch) =>
  postCreateFolder(CONFIG.project_name, parent_folder ? `${parent_folder}/${folder_name}` : folder_name)
    .then(({folders}) => dispatch(addFolders(folders)))
    .catch(error=>alert('Could not create folder!'));

export const removeFolder = ({folder}) => (dispatch) => {
  if (!confirm(`Are you sure you want to remove ${folder.folder_path}?`)) return;

  deleteFolder(
    CONFIG.project_name, folder.folder_path
  )
    .then(({folders}) => dispatch(removeFolders(folders)))
    .catch(error=>alert('Could not remove folder!'));
}

export const protectFolder = ({folder}) => (dispatch) => {
  postProtectFolder(
    CONFIG.project_name, folder.folder_path
  )
    .then(({folders}) => dispatch(addFolders(folders)))
    .catch(error=>alert('Could not protect folder!'));
}

export const unprotectFolder = ({folder}) => (dispatch) => {
  if (!confirm(`Are you sure you want to unprotect ${folder.folder_path}?`)) return;

  postUnprotectFolder(
    CONFIG.project_name, folder.folder_path
  )
    .then(({folders}) => dispatch(addFolders(folders)))
    .catch(error=>alert('Could not unprotect folder!'));
}

export const renameFolder = ({folder, new_folder_path}) => (dispatch) => {
  postRenameFolder(
    CONFIG.project_name, folder.folder_path, new_folder_path
  )
    .then(({folders}) => {
      dispatch(removeFolders([folder]));
      dispatch(addFolders(folders));
    })
    .catch(error=>alert('Could not rename folder!'));
}
