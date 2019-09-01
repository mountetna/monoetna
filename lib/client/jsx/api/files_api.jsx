import { headers, parseJSON, checkStatus,
  json_get, json_delete, json_post } from '../utils/fetch';

export const postRetrieveFiles = (project_name, bucket_name, folder_name) =>
  json_get(`/${project_name}/list/${bucket_name}/${folder_name}`);

export const deleteFile = (project_name, bucket_name, file_name) =>
  json_delete(`/${project_name}/file/remove/${bucket_name}/${file_name}`);

export const postProtectFile = (project_name, bucket_name, file_name) =>
  json_post(`/${project_name}/file/protect/${bucket_name}/${file_name}`);

export const postUnprotectFile = (project_name, bucket_name, file_name) =>
  json_post(`/${project_name}/file/unprotect/${bucket_name}/${file_name}`);

export const postRenameFile = (project_name, bucket_name, file_name, new_file_path) =>
  json_post(`/${project_name}/file/rename/${bucket_name}/${file_name}`, {new_file_path});
