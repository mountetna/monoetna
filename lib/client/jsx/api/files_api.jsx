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
  json_delete(`/${project_name}/remove_file/files/${file_name}`);

export const postProtectFile = (project_name, file_name) =>
  json_post(`/${project_name}/protect_file/files/${file_name}`);

export const postUnprotectFile = (project_name, file_name) =>
  json_post(`/${project_name}/unprotect_file/files/${file_name}`);

export const postRenameFile = (project_name, file_name, new_file_path) =>
  json_post(`/${project_name}/rename_file/files/${file_name}`, {new_file_path});

export const postCreateFolder = (project_name, folder_name) =>
  json_post(`/${project_name}/create_folder/files/${folder_name}`);

export const postProtectFolder = (project_name, folder_name) =>
  json_post(`/${project_name}/protect_folder/files/${folder_name}`);

export const postUnprotectFolder = (project_name, folder_name) =>
  json_post(`/${project_name}/unprotect_folder/files/${folder_name}`);

export const postRenameFolder = (project_name, folder_name, new_folder_path) =>
  json_post(`/${project_name}/rename_folder/files/${folder_name}`, {new_folder_path});

export const deleteFolder = (project_name, folder_name) =>
  json_delete(`/${project_name}/remove_folder/files/${folder_name}`);
