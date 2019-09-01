import { fileKey } from '../utils/file';

// helper to create a new upload entry
const upload = (file, file_name, url) => ({
  file,
  file_name,
  url,
  project_name: CONFIG.project_name,
  file_size: file.size,
  current_byte_position: 0,
  status: 'authorized'
})

const uploads = (old_uploads, action) => {
  if (!old_uploads) old_uploads = {}

  switch(action.type) {
    case 'UPLOAD_SPEED': {
      let { upload, upload_speed } = action;
      let { project_name, file_name } = upload;

      let key = fileKey(upload);

      return {
        ...old_uploads,
        [key]: {
          ...old_uploads[key],
          upload_speed
        }
      };
    }
    case 'UPLOAD_STATUS': {
      let { upload, status } = action;
      let { project_name, file_name, current_byte_position, next_blob_size, next_blob_hash } = upload;
      let key = fileKey(upload);
      return {
        ...old_uploads,
        [key]: {
          ...old_uploads[key],
          ...(status && { status }),
          current_byte_position,
          next_blob_size,
          next_blob_hash
        }
      };
    };
    case 'ADD_UPLOAD': {
      // Copy the selected file data to 'uploads' object.
      let { file, file_name, url } = action;
      let new_upload = upload(file, file_name, url);
      let key = fileKey(new_upload);
      return {
        ...old_uploads,
        [key]: new_upload
      };
    }
    case 'REMOVE_UPLOAD': {
      let { upload } = action;
      let key = fileKey(upload);
      let { [key]: del_upload, ...new_uploads } = old_uploads;

      return new_uploads;
      break;
    }
    default:
      return old_uploads;
  }
};

export default uploads;
