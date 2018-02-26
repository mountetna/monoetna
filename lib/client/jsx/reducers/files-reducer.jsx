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
      return {
        ...state,
        uploads: uploads(state.uploads,action)
      };

    case 'FILE_UPLOAD_RECOVERED':
      setResponseAndIndex(action, fileFails);
      if(index == null) break;

      // Append all of the request items to the local file object.
      fileFails[index] = Object.assign(fileFails[index], response);
      fileFails[index]=Object.assign(action.uploadFile,fileFails[index]);
      uploads.push(fileFails[index]);
      fileFails.splice(index, 1);
      break;

    case 'FILE_INITIALIZED':
    case 'FILE_UPLOAD_ACTIVE':
    case 'FILE_UPLOAD_PAUSED':
      setResponseAndIndex(action, uploads);
      if(index == null) break;

      // Append all of the request items to the local file object.
      uploads[index] = Object.assign(uploads[index], response);
      break;

    case 'FILE_UPLOAD_COMPLETE':
      setResponseAndIndex(action, uploads);
      if(index == null) break;

      /*
       * Move the completed upload metadata from the 'uploads' array to the 
       * 'list' array.
       */
      downloads.push(Object.assign(uploads[index], response));
      uploads.splice(index, 1);
      break;

    case 'FILE_UPLOAD_CANCELLED':
      setResponseAndIndex(action, uploads);
      if(index == null) break;

      // Remove the cancelled upload.
      uploads.splice(index, 1);
      break;

    case 'FILE_REMOVED':
      response = camelCaseIt(action.response.request);

      // Remove the deleted item from uploads.
      index = getMatchingUploadIndex(uploads, response);
      if(index != null) uploads.splice(index, 1);

      // Remove the deleted item from downloads.
      index = getMatchingUploadIndex(downloads, response);
      if(index != null) downloads.splice(index, 1);

      // Remove the deleted item from fileFails.
      index = getMatchingUploadIndex(fileFails, response);
      if(index != null) fileFails.splice(index, 1);
      break;

    case 'CLEAR_UPLOAD':
      for(var a = 0; a < uploads.length; ++a){
        if(action.fileMetadata.reactKey == uploads[a].reactKey){
          uploads.splice(a, 1);
          break;
        }
      }
      break;

    case 'FILE_METADATA_RECEIVED': {
      let { file_list } = action;
      file_list.forEach( file => {
        file = camelCaseIt(file);
        file.reactKey = GENERATE_RAND_KEY();

        if(!file.hasOwnProperty('finishUpload')){
          files.fails.push(file);
        }
        else{
          files.downloads.push(file);
        }
      })
      break;
   }

    case 'QUEUE_UPLOAD':
      for(var a = 0; a < uploads.length; ++a){
        // Remove any 'queued' status
        if(uploads[a].status == 'queued'){
          if(uploads[a].current_byte_position == 0){
            uploads[a].status = 'initialized';
          }
          else{
            uploads[a].status = 'paused';
          }
        }

        // Apply a 'queued' status to a matching file upload.
        if(uploads[a].reactKey == action.reactKey){
          uploads[a].status = 'queued';
        }
      }
      break;

    default:
      return state;
      break;
  }
};

export default files;

