import * as React from 'react';

import ListEntry from './list-entry';
import ListUpload from './list-upload';

export default class ListBody extends React.Component{

  constructor(){

    super();
  }
  
  render(){

    var fileUploads = this['props']['fileUploads'];
    var fileList = this['props']['fileList'];
    
    return (

      <tbody id='list-body-group'>
        
        { (fileUploads.length) ?

            fileUploads.map((fileUpload)=>{

              var redisIndex = fileUpload['redisIndex'];
              return <ListUpload key={ redisIndex } fileUpload={ fileUpload } />
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