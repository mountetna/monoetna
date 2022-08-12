import {
  postRetrieveFiles, postProtectFile, postUnprotectFile, postRenameFile, deleteFile,
  getTouchFile
} from '../api/files_api';
import {errorMessage} from './message_actions';
import {assertIsSome} from "etna-js/utils/asserts";
import DownZip from 'downzip/src/downzip';
import {selectFilesInCurrentFolder} from '../selectors/folder-selector';

const downZip = new DownZip({ scope: document.location.pathname });
const addFiles = (files) => ({type: 'ADD_FILES', files});
const addFolders = (folders) => ({type: 'ADD_FOLDERS', folders});
const removeFiles = (files) => ({type: 'REMOVE_FILES', files});

export const retrieveFiles = ({folder_name, bucket_name}) => (dispatch) =>
  postRetrieveFiles(CONFIG.project_name, bucket_name || 'files', folder_name == undefined ? '' : folder_name)
    .then(({files, folders}) => {
      dispatch(addFiles(files));
      dispatch(addFolders(folders));
    })
    .catch(error => dispatch({type: 'INVALID_FOLDER'}));

export const listFilesRecursive = ({folder_name = "", bucket_name}) => (dispatch) => {
  assertIsSome({folder_name, bucket_name});
  const joinableFolderPath = (folder_name ? folder_name + "/" : "");
  return postRetrieveFiles(CONFIG.project_name, bucket_name, folder_name).then(({files, folders}) =>
    Promise.all(folders.map(({folder_name: subFolder}) => listFilesRecursive({
      folder_name: joinableFolderPath + subFolder,
      bucket_name
    })(dispatch))).then(otherLists => [].concat(files, ...otherLists).sort((a, b) => a.file_path < b.file_path ? -1 : 1))
  );
}

export const downloadFilesZip = ({ folder_name = "", files = [], bucket_name = "" }) => (dispatch) => {
  // a stable identifier to share to cache some downzip work when re-clicked.
  const downloadId = Array.from(files.map(({ file_hash }) => file_hash).join(".")).reduce((s, c) => Math.imul(31, s) + c.charCodeAt(0) | 0, 0)
  let folderParts = folder_name.split("/");
  if (!folder_name) folderParts = [bucket_name];
  const zipName = folderParts.join("--");

  return downZip.downzip(downloadId, zipName, files.map(({ size, download_url, file_path }) => ({
    name: file_path,
    downloadUrl: download_url,
    size
  }))).then(downloadUrl => {
    if (!downloadUrl)  {
      console.error('Could not complete download, your browser may not support Service Workers or the service worker could not be installed.');
      return;
    }

      console.log({ downloadUrl });

      const a = document.createElement('a');
      a.setAttribute('href', downloadUrl);
      document.body.append(a);
      a.click();

      setTimeout(function() {
        a.remove();
      }, 0);
  })
}

export function prepareDownload(files) {
  console.log({files});
}

export const removeFile = ({file, bucket_name}) => (dispatch) => {
  if (!confirm(`Are you sure you want to remove ${file.file_path}?`)) return;

  deleteFile(
    CONFIG.project_name, bucket_name, file.file_path
  )
    .then(({files}) => dispatch(removeFiles(files)))
    .catch(
      errorMessage(dispatch, 'warning', 'File removal failed', error => error)
    );
}

export const protectFile = ({file, bucket_name}) => (dispatch) => {
  postProtectFile(
    CONFIG.project_name, bucket_name, file.file_path
  )
    .then(({files}) => dispatch(addFiles(files)))
    .catch(
      errorMessage(dispatch, 'warning', 'File protection failed', error => error)
    );
}

export const unprotectFile = ({file, bucket_name}) => (dispatch) => {
  if (!confirm(`Are you sure you want to unprotect ${file.file_path}?`)) return;

  postUnprotectFile(
    CONFIG.project_name, bucket_name, file.file_path
  )
    .then(({files}) => dispatch(addFiles(files)))
    .catch(
      errorMessage(dispatch, 'warning', 'File unprotection failed', error => error)
    );
}

export const renameFile = ({file, new_file_path, bucket_name, current_folder, new_bucket_name}) => (dispatch) => {
  postRenameFile(
    CONFIG.project_name, bucket_name, file.file_path, new_file_path, new_bucket_name
  )
    .then(({files}) => {
      dispatch(removeFiles([file]));
      dispatch(addFiles(selectFilesInCurrentFolder({files, current_folder, bucket_name})));
    })
    .catch(
      errorMessage(dispatch, 'warning', 'File renaming failed', error => error)
    );
}

export const touchFile = ({bucket_name, file}) => (dispatch) => {
  getTouchFile(
    CONFIG.project_name, bucket_name, file.file_path
  )
    .then(({files}) => {
      dispatch(removeFiles([file]));
      dispatch(addFiles(files));
    })
    .catch(
      errorMessage(dispatch, 'warning', 'File touching failed', error => error)
    );
}
