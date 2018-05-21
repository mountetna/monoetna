import { headers, parseJSON, checkStatus } from './fetch_utils';

export const postRetrieveFiles = (project_name, folder_name) => {
  return fetch(`/${project_name}/list/files/${folder_name}`, 
  {
    method: 'GET',
    credentials: 'same-origin',
    headers: headers('json')
  }).then(checkStatus).then(parseJSON);
}

export const postCreateFolder = (project_name, folder_name) => {
  return fetch(`/${project_name}/create_folder/files/${folder_name}`,
  {
    method: 'POST',
    credentials: 'same-origin',
    headers: headers('json')
  }).then(checkStatus).then(parseJSON);
}
