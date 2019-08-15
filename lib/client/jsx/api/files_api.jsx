import { headers, parseJSON, checkStatus } from '../utils/fetch';

const json_fetch = (method) => (path, params) => fetch(path,
  {
    method,
    credentials: 'same-origin',
    headers: headers('json'),
    ...params && { body: JSON.stringify(params) }
  }).then(checkStatus).then(parseJSON);

const json_get = json_fetch('GET');
const json_delete = json_fetch('DELETE');
const json_post = json_fetch('POST');

export const postRetrieveFiles = (project_name, bucket_name, folder_name) =>
  json_get(`/${project_name}/list/${bucket_name}/${folder_name}`);

export const postRetrieveBuckets = (project_name) =>
  json_get(`/${project_name}/list`);

export const deleteFile = (project_name, file_name) =>
  json_delete(`/${project_name}/file/remove/files/${file_name}`);

export const postProtectFile = (project_name, file_name) =>
  json_post(`/${project_name}/file/protect/files/${file_name}`);

export const postUnprotectFile = (project_name, file_name) =>
  json_post(`/${project_name}/file/unprotect/files/${file_name}`);

export const postRenameFile = (project_name, file_name, new_file_path) =>
  json_post(`/${project_name}/file/rename/files/${file_name}`, {new_file_path});

export const postCreateFolder = (project_name, folder_name) =>
  json_post(`/${project_name}/folder/create/files/${folder_name}`);

export const postProtectFolder = (project_name, folder_name) =>
  json_post(`/${project_name}/folder/protect/files/${folder_name}`);

export const postUnprotectFolder = (project_name, folder_name) =>
  json_post(`/${project_name}/folder/unprotect/files/${folder_name}`);

export const postRenameFolder = (project_name, folder_name, new_folder_path) =>
  json_post(`/${project_name}/folder/rename/files/${folder_name}`, {new_folder_path});

export const deleteFolder = (project_name, folder_name) =>
  json_delete(`/${project_name}/folder/remove/files/${folder_name}`);
