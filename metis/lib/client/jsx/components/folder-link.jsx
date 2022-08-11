import * as React from 'react';
import { filePath } from 'etna-js/utils/file';

const browse_path = (bucket_name, folder_name, current_folder, log) => {
  let path = filePath(current_folder, folder_name);
  return `/${CONFIG.project_name}/browse/${bucket_name}${path ? `/${path}` : ''}`;
};

export const FolderLink = ({folder_name, log, bucket_name, folder_path}) =>
  <a href={
    browse_path(bucket_name, folder_name, folder_path, log)
  }>{folder_name || bucket_name}</a>;

export const RootFolderLink = ({name}) =>
  <a href={
    `/${CONFIG.project_name}`
  }>{name || CONFIG.project_name}</a>;
