import * as React from 'react';

import ListEntry from './list-entry';
import ListUpload from './list-upload';

export default class ListBody extends React.Component{

  constructor(){

    super();
  }

  startFileUpload(fileUpload){

    this['props'].authorizeFile(fileUpload);
  }
  
  render(){

    var appState = this['props']['appState'];
    var fileUploads = appState['fileUploads'];
    var fileList = appState['fileList'];
    var permissions = appState['userInfo']['permissions'];

    return (

      <tbody id='list-body-group'>
        
        { (fileUploads.length) ?

            fileUploads.map((fileUpload)=>{

              var listUpload = {

                'key': fileUpload['redisIndex'],
                'fileUpload': fileUpload,
                'permissions': permissions,
                'callbacks': {

                  'startFileUpload': this['startFileUpload'].bind(this)
                }
              }
              return <ListUpload { ...listUpload } />
            })
          : '' }

        { (fileList.length) ?
            
            fileList.map((fileInfo)=>{

              var redisIndex = fileInfo['redisIndex'];
              return <ListEntry key={ redisIndex } fileInfo={ fileInfo } />
            })
          : '' }
      </tbody>
    );
  }
}