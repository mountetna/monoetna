// helper to create a new file entry
const upload = (file, file_name, url) => ({
  file,
  file_name,
  url,
  project_name: CONFIG.project_name,
  file_size: file.size,
  current_byte_position: 0,
  status: 'authorized'
})

const file_key = ({project_name, file_name}) => `${project_name}:${file_name}`;

const uploads = (old_uploads, action) => {
  if (!old_uploads) old_uploads = {}

  switch(action.type) {
    case 'FILE_UPLOAD_SPEED': {
      let { upload, upload_speed } = action;
      let { project_name, file_name } = upload;

      let key = file_key(upload);

      return {
        ...old_uploads,
        [key]: {
          ...old_uploads[key],
          upload_speed
        }
      };
    }
    case 'FILE_UPLOAD_STATUS': {
      let { upload, status } = action;
      let { project_name, file_name, current_byte_position, next_blob_size, next_blob_hash } = upload;
      let key = file_key(upload);
      return {
        ...old_uploads,
        [key]: {
          ...old_uploads[key],
          status,
          current_byte_position,
          next_blob_size,
          next_blob_hash
        }
      };
    };
    case 'FILE_UPLOAD_SELECT_PROJECT': {
      let { upload, permission: { project_name, role } } = action;
      let key = file_key(upload);
      return {
        ...old_uploads,
        [key]: {
          ...old_uploads[key],
          project_name, role
        }
      };
    };
    case 'FILE_UPLOAD_AUTHORIZED': {
      // Copy the selected file data to 'uploads' object.
      let { file, file_name, url } = action;
      let new_upload = upload(file, file_name, url);
      let key = file_key(new_upload);
      return {
        ...old_uploads,
        [key]: new_upload
      };
    }
    case 'FILE_UPLOAD_REMOVED': {
      let { upload, status } = action;
      let { project_name, file_name, current_byte_position, next_blob_size, next_blob_hash } = upload;
      let key = file_key(upload);
      let { [key]: del_upload, ...new_uploads } = old_uploads;

      return new_uploads;
      break;
    }
    default:
      return old_uploads;
  }
};

const folders = (old_folders, action) => {
  if (!old_folders) old_folders = {};

  switch(action.type) {
    case 'ADD_FOLDERS': {
      let { folders } = action;

      let new_folders = folders.reduce((c, file) => {
        let key = file_key(file);
        c[key] = file;
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

const files = (old_files, action) => {
  if (!old_files) old_files = {};

  switch(action.type) {
    case 'ADD_FILES': {
      let { files } = action;

      let new_files = files.reduce((c, file) => {
        let key = file_key(file);
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

const directory = (state, action) => {
  if (!state) state = {
    folders: {},
    files: {},
    uploads: {},
    fails: [],
    current_folder: null
  };

  switch(action.type) {
    case 'FILE_UPLOAD_STATUS':
    case 'FILE_UPLOAD_SPEED':
    case 'FILE_UPLOAD_AUTHORIZED':
    case 'FILE_UPLOAD_SELECT_PROJECT':
    case 'FILE_UPLOAD_REMOVED':
      return {
        ...state,
        uploads: uploads(state.uploads,action)
      };
    case 'ADD_FOLDERS':
      return {
        ...state,
        folders: folders(state.folders,action)
      };

    case 'ADD_FILES':
      return {
        ...state,
        files: files(state.files,action)
      };

    case 'SET_CURRENT_FOLDER':
      return {
        ...state,
        current_folder: action.folder_name
      }
    default:
      return state;
      break;
  }
};

export default directory;

