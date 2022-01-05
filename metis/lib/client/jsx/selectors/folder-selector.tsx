import {Folder} from '../types/metis-types';

export const selectFoldersInCurrentFolder = ({
  folders,
  current_folder,
  bucket_name
}: {
  folders: Folder[];
  current_folder: string;
  bucket_name: string;
}) => {
  return folders.filter((folder) => {
    if (bucket_name && folder.bucket_name !== bucket_name) return false;
    if (current_folder && !folder.folder_path.startsWith(current_folder))
      return false;
    if (current_folder === '' && folder.folder_path.includes('/')) return false;

    return true;
  });
};
