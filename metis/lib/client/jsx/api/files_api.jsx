import { json_get, json_delete, json_post } from 'etna-js/utils/fetch';

// TODO: These really ought to be refactored into a service class that knows about metis' domain
// in various environments, rather than assuming the location of the currently viewed page contains
// these endpoints.
export const postRetrieveFiles = (project_name, bucket_name, folder_name) =>
  json_get(`${window.location.origin}/${project_name}/list/${bucket_name}/${folder_name}`);

export const deleteFile = (project_name, bucket_name, file_name) =>
  json_delete(`${window.location.origin}/${project_name}/file/remove/${bucket_name}/${file_name}`);

export const postProtectFile = (project_name, bucket_name, file_name) =>
  json_post(`${window.location.origin}/${project_name}/file/protect/${bucket_name}/${file_name}`);

export const postUnprotectFile = (project_name, bucket_name, file_name) =>
  json_post(`${window.location.origin}/${project_name}/file/unprotect/${bucket_name}/${file_name}`);

export const postRenameFile = (project_name, bucket_name, file_name, new_file_path, new_bucket_name) => {
  let payload = {
    new_file_path
  };

  if (new_bucket_name) payload.new_bucket_name = new_bucket_name;
  return json_post(`${window.location.origin}/${project_name}/file/rename/${bucket_name}/${file_name}`, payload);
}


export const getTouchFile = (project_name, bucket_name, file_name) =>
  json_get(`${window.location.origin}/${project_name}/file/touch/${bucket_name}/${file_name}`);
