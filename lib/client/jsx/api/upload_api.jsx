import { headers, checkStatus } from './fetch_utils';

export const postAuthorizeUpload = (upload) => {
  let request = {
    project_name: upload.projectName,
    file_name: upload.fileName
  };

  return fetch('/authorize/upload', 
  {
    method: 'POST',
    credentials: 'same-origin',
    headers: headers('json'),
    body: JSON.stringify(request)
  }).then(checkStatus);
}
