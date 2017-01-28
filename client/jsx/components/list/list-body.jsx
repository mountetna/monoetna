import * as React from 'react';

import ListEntry from './list-entry';
import ListUpload from './list-upload';
import ListUploadFailed from './list-upload-failed';

export default class ListBody extends React.Component{

  constructor(){

    super();
  }
  
  render(){

    var fileUploads = this['props']['fileData']['fileUploads'];
    var fileList = this['props']['fileData']['fileList'];
    var fileFails = this['props']['fileData']['fileFails'];
    var permissions = this['props']['userInfo']['permissions'];

    return (

      <tbody id='list-body-group'>

        {/* Render the failed uploads. */}
        { (fileFails['length']) ? 

            fileFails.map((failedFile)=>{

              var failedUpload = {

                'key': failedFile['redisIndex'],
                'failedFile': failedFile,
                'callbacks': {

                  'recoverUpload': this['props']['recoverUpload'],
                  'removeFile': this['props']['removeFile']
                }
              };

              return <ListUploadFailed { ...failedUpload } />;
            })
          : '' }

        {/* Render the incomplete uploads. */}
        { (fileUploads['length']) ?

            fileUploads.map((fileUpload)=>{

              var listUpload = {

                'key': fileUpload['redisIndex'],
                'fileUpload': fileUpload,
                'permissions': permissions,
                'callbacks': {

                  'initializeUpload': this['props']['initializeUpload'],
                  'queueUpload': this['props']['queueUpload'],
                  'pauseUpload': this['props']['pauseUpload'],
                  'cancelUpload': this['props']['cancelUpload']
                }
              };

              return <ListUpload { ...listUpload } />
            })
          : '' }

        {/* Render the complete uploads. */}
        { (fileList['length']) ?
            
            fileList.map((fileInfo)=>{

              var listEntry = {

                'key': fileInfo['redisIndex'],
                'fileInfo': fileInfo,
                'callbacks': {

                  'removeFile': this['props']['removeFile']
                }
              };

              return <ListEntry { ...listEntry } />
            })
          : '' }
      </tbody>
    );
  }
}