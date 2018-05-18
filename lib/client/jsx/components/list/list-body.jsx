import * as React from 'react';
import * as ReactRedux from 'react-redux';

import ListEntry from './list-entry';
import ListUpload from './list-upload';
import ListUploadFailed from './list-upload-failed';

class ListBody extends React.Component{
  constructor() {
    super();
  }
  
  render() {
    let { files, widths, user } = this.props;
    let { uploads, downloads, fails } = files;
    let { permissions } = user;

    let { recoverUpload, removeFailed, initializeUpload,
      queueUpload, pauseUpload, cancelUpload, removeFile,
      clearUpload, selectProject } = this.props;

    return (
      <div id='list-body-group'>
        {/* Render the failed uploads. */}
        { (fails.length) ? 
            fails.map((failedFile)=>{
              let failedUpload = {
                key: 'file-failed-'+failedFile.reactKey,
                failedFile: failedFile,
                callbacks: { recoverUpload, removeFailed }
              };

              return <ListUploadFailed { ...failedUpload } />;
            })
          : '' }
        {/* Render the incomplete uploads. */}
        { (Object.keys(uploads).length) ?
            Object.keys(uploads).map((key)=>{
              let upload = uploads[key];
              console.log(upload);
              let listUpload = {
                key,
                upload,
                user,
                permissions,
                callbacks: {
                  initializeUpload, queueUpload, pauseUpload,
                  cancelUpload, removeFile, clearUpload, selectProject
                }
              };

              return <ListUpload { ...listUpload } />
            })
          : '' }
        {/* Render the complete uploads. */}
        { (Object.keys(downloads).length) ?
            Object.keys(downloads).map(key=>{
              let file = downloads[key];
              let listEntry = {
                key,
                file,
                widths,
                callbacks: { removeFile }
              };

              return <ListEntry { ...listEntry } />
            })
          : '' }
      </div>
    );
  }
}


const ListBodyContainer = ReactRedux.connect(
  // map state
  ({user,files}) => ({user,files}),

  // map dispatch
  (dispatch) => ({
    selectProject: (upload, permission) => dispatch({ type: 'FILE_UPLOAD_SELECT_PROJECT', upload, permission }),
    initializeUpload: (upload) => dispatch({ type: 'AUTHORIZE_UPLOAD', upload }),
    queueUpload: (upload) => dispatch({ type: 'QUEUE_UPLOAD', upload }),
    pauseUpload: (upload) => dispatch({ type: 'PAUSE_UPLOAD', upload }),
    cancelUpload: (upload) => dispatch({ type: 'CANCEL_UPLOAD', upload }),
    recoverUpload: (uploadFile, fileMetadata) => dispatch({ type: 'RECOVER_UPLOAD', uploadFile, fileMetadata }),
    removeFile: (fileMetadata) => dispatch({ type: 'REMOVE_FILE', fileMetadata }),
    clearUpload: (fileMetadata) => dispatch({ type: 'CLEAR_UPLOAD', fileMetadata }),
    removeFailed: (fileMetadata) => dispatch({ type: 'REMOVE_FAILED', fileMetadata })
  })
)(ListBody);

export default ListBodyContainer;
