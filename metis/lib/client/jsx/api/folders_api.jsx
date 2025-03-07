import {json_get, json_delete, json_post} from 'etna-js/utils/fetch';

export const postCreateFolder = (project_name, bucket_name, folder_name) =>
  json_post(`/${project_name}/folder/create/${bucket_name}/${folder_name}`);

export const postProtectFolder = (project_name, bucket_name, folder_name) =>
  json_post(`/${project_name}/folder/protect/${bucket_name}/${folder_name}`);

export const postUnprotectFolder = (project_name, bucket_name, folder_name) =>
  json_post(`/${project_name}/folder/unprotect/${bucket_name}/${folder_name}`);

export const postRenameFolder = (
  project_name,
  bucket_name,
  folder_name,
  new_folder_path,
  new_bucket_name = null
) => {
  let payload = {
    new_folder_path
  };

  if (new_bucket_name) payload.new_bucket_name = new_bucket_name;

  return json_post(
    `/${project_name}/folder/rename/${bucket_name}/${folder_name}`,
    payload
  );
};

export const deleteFolder = (
  project_name,
  bucket_name,
  folder_name,
  recursive
) =>
  json_delete(
    `/${project_name}/folder/remove/${bucket_name}/${folder_name}`,
    recursive ? {recursive: true} : {}
  );

export const getTouchFolder = (project_name, bucket_name, folder_name) =>
  json_get(`/${project_name}/folder/touch/${bucket_name}/${folder_name}`);

export const postCopyFolder = (
  project_name,
  bucket_name,
  folder_name,
  new_parent_folder,
  new_bucket_name = null
) => {
  let payload = {
    new_parent_folder
  };

  if (new_bucket_name) payload.new_bucket_name = new_bucket_name;

  return json_post(
    `/${project_name}/folder/copy/${bucket_name}/${folder_name}`,
    payload
  );
};
