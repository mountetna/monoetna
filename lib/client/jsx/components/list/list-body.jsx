import * as React from 'react';

import ListEntry from './list-entry';
import ListUpload from './list-upload';
import ListUploadFailed from './list-upload-failed';

export default class ListBody extends React.Component{
  constructor() {
    super();
  }
  
  render() {
    let { files, user } = this.props;
    let { uploads, downloads, fails } = files;
    let { permissions } = user;

    let { recoverUpload, removeFailed, initializeUpload,
      queueUpload, pauseUpload, cancelUpload, removeFile,
      clearUpload, selectProject } = this.props;

    return (
      <tbody id='list-body-group'>
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
                callbacks: { removeFile }
              };

              return <ListEntry { ...listEntry } />
            })
          : '' }
      </tbody>
    );
  }
}
