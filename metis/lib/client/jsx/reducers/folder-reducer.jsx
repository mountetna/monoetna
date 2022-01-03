import { folderKey } from 'etna-js/utils/file';

const addFolders = (action, old_folders) => {
  let { folders, bucket_name } = action;

  let new_folders = folders.reduce((c, folder) => {
    let key = folderKey(folder);

    if (bucket_name && folder.bucket_name !== bucket_name) return c;
    
    c[key] = folder;
    return c;
  }, {});

  return {
    ...old_folders,
    ...new_folders
  };
}

const removeFolders = (action, old_folders) => {
  let { folders } = action;

  let bad_folders = folders.reduce((c, folder) => {
    let key = folderKey(folder);
    c[key] = true;
    return c;
  }, {});

  let new_folders = Object.keys(old_folders).reduce((c, key) => {
    if (!(key in bad_folders)) c[key] = old_folders[key];
    return c;
  }, {});

  return new_folders;
}

const folders = (old_folders, action) => {
  if (!old_folders) old_folders = {};

  switch(action.type) {
    case 'ADD_FOLDERS':
      return addFolders(action, old_folders);
    case 'REMOVE_FOLDERS':
      return removeFolders(action, old_folders);
    default:
      return old_folders;
  }
};

export default folders;
