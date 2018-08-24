export const fileKey = ({project_name, file_name}) => `${project_name}:${file_name}`;
export const folderKey = ({project_name, folder_name}) => `${project_name}:${folder_name}`;

export const filePath = (current_folder, file_name) =>
  `${
    current_folder == undefined ? '' : `${current_folder}/`
  }${file_name}`;
