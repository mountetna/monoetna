import { fileKey } from 'etna-js/utils/file';

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
};

const addFiles = (old_files, action) => {
  let { files } = action;

  let new_files = files.reduce((c, file) => {
    let key = fileKey(file);
    c[key] = file;
    return c;
  }, {});

  return {
    ...old_files,
    ...new_files
  };
};

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
