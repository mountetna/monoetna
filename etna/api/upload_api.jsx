import { form_post, json_post } from '../utils/fetch';

export const postAuthorizeUpload = (base_url, project_name, bucket_name, file_path) => {
  let request = { project_name, file_path, bucket_name };

  return json_post(base_url + '/authorize/upload', request)
}

export const postUploadStart = (upload_url, request) => json_post(upload_url, { action: 'start', ...request });
export const postUploadCancel = (upload_url, request) => json_post(upload_url, { action: 'cancel', ...request });
export const postUploadBlob = (upload_url, request) => form_post(upload_url, { action: 'blob', ...request });
