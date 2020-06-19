export const fileKey = ({project_name, file_name}) => `${project_name}:${file_name}`;
export const folderKey = ({project_name, folder_name}) => `${project_name}:${folder_name}`;

export const filePath = (folder_name, file_name) =>
  `${
    !folder_name ? '' : `${folder_name}/`
  }${!file_name ? '' : file_name}`;
