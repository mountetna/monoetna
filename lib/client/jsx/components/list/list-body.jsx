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

    return (
      <tbody id='list-body-group'>
        {/* Render the failed uploads. */}
        { (fails.length) ? 
            fails.map((failedFile)=>{
              let failedUpload = {
                key: 'file-failed-'+failedFile.reactKey,
                failedFile: failedFile,
                callbacks: {
                  recoverUpload: this.props.recoverUpload,
                  removeFailed: this.props.removeFailed,
                }
              };

              return <ListUploadFailed { ...failedUpload } />;
            })
          : '' }
        {/* Render the incomplete uploads. */}
        { (Object.keys(uploads).length) ?
            Object.keys(uploads).map((key)=>{
              let upload = uploads[key];
              let listUpload = {
                key: 'file-upload-'+key,
                key,
                upload,
                permissions,
                callbacks: {
                  initializeUpload: this.props.initializeUpload,
                  queueUpload: this.props.queueUpload,
                  pauseUpload: this.props.pauseUpload,
                  cancelUpload: this.props.cancelUpload,
                  removeFile: this.props.removeFile,
                  clearUpload: this.props.clearUpload,
                  selectProject: this.props.selectProject
                }
              };

              return <ListUpload { ...listUpload } />
            })
          : '' }
        {/* Render the complete uploads. */}
        { (downloads.length) ?
            downloads.map((file, i)=>{
              let listEntry = {
                key: i,
                file,
                callbacks: {
                  removeFile: this.props.removeFile
                }
              };

              return <ListEntry { ...listEntry } />
            })
          : '' }
      </tbody>
    );
  }
}
