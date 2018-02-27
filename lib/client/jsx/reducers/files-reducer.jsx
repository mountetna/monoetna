// helper to create a new file entry
const upload = (file, url) => ({
  file,
  file_name: file.name,
  url,
  project_name: CONFIG.project_name,
  file_size: file.size,
  current_byte_position: 0,
  status: 'authorized'
})

const upload_key = ({project_name, file_name}) => `${project_name}:${file_name}`;

const uploads = (old_uploads, action) => {
  if (!old_uploads) old_uploads = {}

  switch(action.type) {
    case 'FILE_UPLOAD_SPEED': {
      let { upload, upload_speed } = action;
      let { project_name, file_name } = upload;

      let key = upload_key(upload);

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
      let key = upload_key(upload);
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
      let key = upload_key(upload);
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
      let { file, url } = action;
      let new_upload = upload(file, url);
      let key = upload_key(new_upload);
      return {
        ...old_uploads,
        [key]: new_upload
      };
    }
    case 'FILE_UPLOAD_REMOVED': {
      let { upload, status } = action;
      let { project_name, file_name, current_byte_position, next_blob_size, next_blob_hash } = upload;
      let key = upload_key(upload);
      let { [key]: del_upload, ...new_uploads } = old_uploads;

      return new_uploads;
      break;
    }
    default:
      return old_uploads;
  }
};

const files = (state, action) => {
  if (!state) state = {
    downloads: [],
    uploads: {},
    fails: []
  } 

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

    default:
      return state;
      break;
  }
};

export default files;

