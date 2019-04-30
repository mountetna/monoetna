import * as React from 'react';

export const FolderLink = ({folder_name, folder_path}) =>
  <a href={
    `/${CONFIG.project_name}/browse/${
      folder_path ? `${folder_path}/` : ''
    }${folder_name}`
  }>{folder_name}</a>;

export const RootFolderLink = ({name}) =>
  <a href={
    `/${CONFIG.project_name}`
  }>{name || CONFIG.project_name}</a>;
