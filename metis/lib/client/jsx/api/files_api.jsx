import {json_get, json_delete, json_post} from 'etna-js/utils/fetch';

export const postRetrieveFiles = (project_name, bucket_name, folder_name) =>
  json_get(`/${project_name}/list/${bucket_name}/${folder_name}`);

export const deleteFile = (project_name, bucket_name, file_name) =>
  json_delete(`/${project_name}/file/remove/${bucket_name}/${file_name}`);

export const postProtectFile = (project_name, bucket_name, file_name) =>
  json_post(`/${project_name}/file/protect/${bucket_name}/${file_name}`);

export const postUnprotectFile = (project_name, bucket_name, file_name) =>
  json_post(`/${project_name}/file/unprotect/${bucket_name}/${file_name}`);

export const postRestrictFile = (project_name, bucket_name, file_name) =>
  json_post(`/${project_name}/file/restrict/${bucket_name}/${file_name}`);

export const postUnrestrictFile = (project_name, bucket_name, file_name) =>
  json_post(`/${project_name}/file/unrestrict/${bucket_name}/${file_name}`);

export const postRenameFile = (
  project_name,
  bucket_name,
  file_name,
  new_file_path,
  new_bucket_name
) => {
  let payload = {
    new_file_path
  };

  if (new_bucket_name) payload.new_bucket_name = new_bucket_name;
  return json_post(
    `/${project_name}/file/rename/${bucket_name}/${file_name}`,
    payload
  );
};

export const getTouchFile = (project_name, bucket_name, file_name) =>
  json_get(`/${project_name}/file/touch/${bucket_name}/${file_name}`);

export const postCopyFiles = (project_name, revisions) => {
  let payload = {
    revisions
  };

  return json_post(`/${project_name}/files/copy`, payload);
};
