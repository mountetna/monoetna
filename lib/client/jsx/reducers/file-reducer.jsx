import { fileKey } from '../utils/file';

const files = (old_files, action) => {
  if (!old_files) old_files = {};

  switch(action.type) {
    case 'ADD_FILES': {
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
    }
    default:
      return old_files;
  }
};

export default files;
