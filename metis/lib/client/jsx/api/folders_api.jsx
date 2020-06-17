import { json_get, json_delete, json_post } from '../utils/fetch';

export const postCreateFolder = (project_name, bucket_name, folder_name) =>
  json_post(`/${project_name}/folder/create/${bucket_name}/${folder_name}`);

export const postProtectFolder = (project_name, bucket_name, folder_name) =>
  json_post(`/${project_name}/folder/protect/${bucket_name}/${folder_name}`);

export const postUnprotectFolder = (project_name, bucket_name, folder_name) =>
  json_post(`/${project_name}/folder/unprotect/${bucket_name}/${folder_name}`);

export const postRenameFolder = (project_name, bucket_name, folder_name, new_folder_path) =>
  json_post(`/${project_name}/folder/rename/${bucket_name}/${folder_name}`, {new_folder_path});

export const deleteFolder = (project_name, bucket_name, folder_name) =>
  json_delete(`/${project_name}/folder/remove/${bucket_name}/${folder_name}`);
