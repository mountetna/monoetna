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

export const deleteFolder = (project_name, folder_name) => {
  return fetch(`/${project_name}/remove_folder/files/${folder_name}`,
  {
    method: 'DELETE',
    credentials: 'same-origin',
    headers: headers('json')
  }).then(checkStatus).then(parseJSON);
}

export const deleteFile = (project_name, file_name) => {
  return fetch(`/${project_name}/remove_file/files/${file_name}`,
  {
    method: 'DELETE',
    credentials: 'same-origin',
    headers: headers('json')
  }).then(checkStatus).then(parseJSON);
}

export const postProtectFile = (project_name, file_name) => {
  return fetch(`/${project_name}/protect_file/files/${file_name}`,
  {
    method: 'POST',
    credentials: 'same-origin',
    headers: headers('json')
  }).then(checkStatus).then(parseJSON);
}

export const postUnprotectFile = (project_name, file_name) => {
  return fetch(`/${project_name}/unprotect_file/files/${file_name}`,
  {
    method: 'POST',
    credentials: 'same-origin',
    headers: headers('json')
  }).then(checkStatus).then(parseJSON);
}
