import { folderKey } from '../utils/file';

const folders = (old_folders, action) => {
  if (!old_folders) old_folders = {};

  switch(action.type) {
    case 'ADD_FOLDERS': {
      let { folders } = action;

      let new_folders = folders.reduce((c, folder) => {
        let key = folderKey(folder);
        c[key] = folder;
        return c;
      }, {});

      return {
        ...old_folders,
        ...new_folders
      };
    }
    default:
      return old_folders;
  }
};

export default folders;
