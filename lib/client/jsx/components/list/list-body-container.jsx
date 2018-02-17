import * as ReactRedux from 'react-redux';
import ListBody from './list-body';

const mapStateToProps = (state, ownProps)=>{
  return {
    user: state.user,
    files: state.files
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{
  return {
    selectProject: (upload, permission) => {
      var action = { type: 'FILE_UPLOAD_SELECT_PROJECT', upload, permission };
      dispatch(action);
    },

    initializeUpload: (upload)=>{
      var action = { type: 'AUTHORIZE_UPLOAD', upload };
      dispatch(action);
    },

    queueUpload: (upload)=>{
      var action = { type: 'QUEUE_UPLOAD', upload };
      dispatch(action);
    },

    pauseUpload: ()=>{
      var action = { type: 'PAUSE_UPLOAD' };
      dispatch(action);
    },

    cancelUpload: ()=>{
      var action = { type: 'CANCEL_UPLOAD' };
      dispatch(action);
    },

    recoverUpload: (uploadFile, fileMetadata)=>{
      var action = { 

        type: 'RECOVER_UPLOAD',
        uploadFile,
        fileMetadata
      };

      dispatch(action);
    },

    removeFile: (fileMetadata)=>{
      var action = { type: 'REMOVE_FILE', fileMetadata };
      dispatch(action);
    },

    /* 
     * This function removes a file/upload from the local data. When a file has
     * not yet been authorized on the server we just remove it from the local
     * store.
     */
    clearUpload: (fileMetadata)=>{
      var action = { type: 'CLEAR_UPLOAD', fileMetadata };
      dispatch(action);
    },

    removeFailed: (fileMetadata)=>{
      var action = { type: 'REMOVE_FAILED', fileMetadata };
      dispatch(action);
    }
  };
}

const ListBodyContainer = ReactRedux.connect(
  mapStateToProps,
  mapDispatchToProps,
)(ListBody);

export default ListBodyContainer;
