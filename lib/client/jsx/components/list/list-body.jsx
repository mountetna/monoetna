import * as React from 'react';
import * as ReactRedux from 'react-redux';

import { ListEntry, ListFolder } from './list-entry';
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

    let download_files = Object.values(downloads).filter(f=>!f.is_folder);
    let download_folders = Object.values(downloads).filter(f=>f.is_folder);

    let { recoverUpload, removeFailed, initializeUpload,
      queueUpload, pauseUpload, cancelUpload, removeFile,
      clearUpload, selectProject } = this.props;

    return (
      <div id='list-body-group'>
        {/* Render the incomplete uploads. */}
        { (Object.keys(uploads).length) ?
            Object.keys(uploads).map((key)=>{
              let upload = uploads[key];
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
        {
          (download_folders.length) ?
            download_folders.map(folder=>{
              let listEntry = {
                key: folder.file_name,
                file: folder,
                widths,
                callbacks: { removeFile }
              };

              return <ListFolder { ...listEntry } />
            })
            : '' 
        }
        {
          (download_files.length) ?
            download_files.map(file=>{
              let listEntry = {
                key: file.file_name,
                file,
                widths,
                callbacks: { removeFile }
              };

              return <ListEntry { ...listEntry } />
            })
            : '' 
        }
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
