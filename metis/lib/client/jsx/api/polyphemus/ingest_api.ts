import {json_get, json_post} from 'etna-js/utils/fetch';

export const listDirectory = (host: string, directory: string) =>
  json_get(
    `${CONFIG.polyphemus_host}/${CONFIG.project_name}/ingest/list/${host}/${directory}`
  );

export const enqueueDirectoryFiles = (host: string, directory: string) =>
  json_post(
    `${CONFIG.polyphemus_host}/${CONFIG.project_name}/ingest/enqueue/${host}/${directory}`
  );
