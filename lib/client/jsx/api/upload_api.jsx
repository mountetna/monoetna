import { headers, parseJSON, checkStatus } from './fetch_utils';

export const postAuthorizeUpload = ({ project_name, file_name }) => {
  let request = { project_name, file_name };

  return fetch('/authorize/upload', 
  {
    method: 'POST',
    credentials: 'same-origin',
    headers: headers('json'),
    body: JSON.stringify(request)
  }).then(checkStatus);
}

const postUpload = (upload_url, request) => {
  return fetch(upload_url,
  {
    method: 'POST',
    credentials: 'same-origin',
    headers: headers('json'),
    body: JSON.stringify(request)
  })
}

export const postUploadStart = (upload_url, request) => {
  return postUpload(upload_url, { ...request, action: 'start' }).then(checkStatus).then(parseJSON);
}

export const postUploadBlob = ({ upload, blob, new_blob_size, new_blob_hash}) => {
  let form = new FormData();
  let { url } = upload;

  form.append('action', 'blob');
  form.append('blob_data', blob);
  form.append('next_blob_size', new_blob_size);
  form.append('next_blob_hash', new_blob_hash);

  return fetch(url,
  {
    method: 'POST',
    credentials: 'same-origin',
    body: form
  }).then(checkStatus).then(parseJSON);
}
