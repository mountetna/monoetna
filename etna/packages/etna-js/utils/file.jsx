import {assertIsSome} from "./asserts";

export const fileKey = ({project_name, file_name}) => `${assertIsSome({project_name})}:${assertIsSome({file_name})}`;

export const folderKey = ({project_name, folder_name}) => `${assertIsSome({project_name})}:${assertIsSome({folder_name})}`;

export const filePath = (folder_name, file_name) =>
  `${
    !folder_name ? '' : `${folder_name}/`
  }${!file_name ? '' : file_name}`;

export const includesFolders = (path) => path && path.includes("/");