import {json_get, json_post} from 'etna-js/utils/fetch';

export const listDirectory = (host: string, directory: string) =>
  json_get(
    `${CONFIG.polyphemus_host}/${CONFIG.project_name}/list/${host}/${directory}`
  );

export const queueDirectoryFiles = (host: string, directory: string) =>
  json_post(
    `${CONFIG.polyphemus_host}/${CONFIG.project_name}/queue/${host}/${directory}`
  );
