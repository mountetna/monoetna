import { headers, checkStatus } from './fetch_utils';

export const uploadAuthorize = (upload) => {
  let request = {
    project_name: upload.projectName,
    file_name: upload.fileName
  };

  return fetch('/authorize/upload', 
  {
    method: 'POST',
    credentials: 'same-origin',
    headers: headers('json'),
    sendType: 'serial',
    body: JSON.stringify(request)
  }).then(checkStatus);
}
