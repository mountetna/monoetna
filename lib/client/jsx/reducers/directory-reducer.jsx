import uploads from './upload-reducer';
import folders from './folder-reducer';
import files from './file-reducer';

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

    case 'INVALID_FOLDER':
      return {
        ...state,
        current_folder: '\ninvalid\n'
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

