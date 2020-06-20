import { fileKey } from '../utils/file';

const removeFiles = (old_files, action) => {
  let { files } = action;

  let removed_keys = files.reduce(
    (rk,file) => {
      rk[fileKey(file)] = true;
      return rk;
    }, {}
  );

  return Object.keys(old_files).filter(key => !(key in removed_keys)).reduce(
    (nf,key) => { nf[key] = old_files[key]; return nf; }, {}
  );
}

const addFiles = (old_files, action) => {
  let { files } = action;

  let new_files = files.reduce((c, file) => {
    // Recursive files won't end up inside the current directory view.
    // TODO: Consolidate the file reducer logic here with etna-js so that more
    // logic can be shared and a smarter interface designed.
    if (!file.webkitRelativePath || file.webkitRelativePath == file.name) {
      let key = fileKey(file);
      c[key] = file;
    }
    return c;
  }, {});

  return {
    ...old_files,
    ...new_files
  };
}

const files = (old_files, action) => {
  if (!old_files) old_files = {};

  switch(action.type) {
    case 'REMOVE_FILES':
      return removeFiles(old_files, action);
    case 'ADD_FILES':
      return addFiles(old_files,action);
    default:
      return old_files;
  }
};

export default files;
