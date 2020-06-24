import uploads from 'etna-js/upload/reducers/upload-reducer';
import folders from './folder-reducer';
import files from './file-reducer';
import buckets from './bucket-reducer';

import {
  UPLOAD_STATUS,
  UPLOAD_SPEED,
  ADD_UPLOAD,
  REMOVE_UPLOAD
} from 'etna-js/upload/actions/upload_actions';

const directory = (state, action) => {
  if (!state) state = {
    folders: {},
    files: {},
    buckets: [],
    uploads: {},
    fails: [],
    current_folder: null
  };

  switch (action.type) {
    case UPLOAD_STATUS:
    case UPLOAD_SPEED:
    case ADD_UPLOAD:
    case REMOVE_UPLOAD:
      return {
        ...state,
        uploads: uploads(state.uploads, action)
      };
    case 'ADD_FOLDERS':
    case 'REMOVE_FOLDERS':
      return {
        ...state,
        folders: folders(state.folders, action)
      };

    case 'ADD_BUCKETS':
    case 'REMOVE_BUCKET':
      return {
        ...state,
        buckets: buckets(state.buckets, action)
      };

    case 'ADD_FILES':
    case 'REMOVE_FILES':
      return {
        ...state,
        files: files(state.files, action)
      };

    case 'INVALID_FOLDER':
      return {
        ...state,
        current_folder: '\ninvalid\n'
      };
    case 'SET_CURRENT_FOLDER':
      return {
        ...state,
        current_folder: action.folder_name,
        current_bucket: action.bucket_name
      }
    default:
      return state;
      break;
  }
};

export default directory;

