import { headers, parseJSON, checkStatus } from './fetch_utils';

export const postRetrieveFiles = (project_name) => {
  return fetch(`/${project_name}/files`, 
  {
    method: 'GET',
    credentials: 'same-origin',
    headers: headers('json')
  }).then(checkStatus).then(parseJSON);
}

export const postCreateFolder = (project_name, folder_name) => {
  return fetch(`/${project_name}/folder/${folder_name}`, 
  {
    method: 'POST',
    credentials: 'same-origin',
    headers: headers('json')
  }).then(checkStatus).then(parseJSON);
}
