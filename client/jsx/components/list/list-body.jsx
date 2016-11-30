import * as React from 'react';

import ListEntry from './list-entry';
import ListUpload from './list-upload';

export default class ListBody extends React.Component{

  constructor(){

    super();
  }

  updateFileUpload(fileUpload){

    console.log(fileUpload);

    /*
    1. Check the new file name...
    2. If good then set...
    3. If no good then say so and reset.
    */
  }

  startFileUpload(fileUpload){

    this['props'].authorizeFile(fileUpload);
  }
  
  render(){

    var fileUploads = this['props']['fileUploads'];
    var fileList = this['props']['fileList'];
    var permissions = this['props']['permissions'];
    
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