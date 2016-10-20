import * as React from 'react'

import UploadMeter from './upload-meter'
import UploadControl from './upload-control'

export default class ListUpload extends React.Component{

  constructor(){

    super();
  }

  parseFileName(){

    var file = this['props']['fileUpload'];
    var origName = file['original_name'];
    var fileName = file['file_name'];

    if(origName != fileName){ 

      return (

        <span>

          { fileName }
          <span className='light-text'>

            { ' ('+ origName +')' }
          </span>
        </span>
      );
    }

    return fileName;
  }

  parseFileStatus(){

    var file = this['props']['fileUpload'];
    var status = file['status'];
    var date = PARSE_TIMESTAMP(Date.now() / 1000);
    var user = file['user_email'];
    switch(file['status']){

      case 'queued':

        status = 'File upload queued.' 
        break;
      case 'initialized':

        status = 'File upload initialized by '+ user;
        break;
      case 'authorized':

        status = 'File upload authorized on '+ date;
        break;
      case 'active':

        status = 'File uploading...'
        break;
      case 'complete':

        status = 'Uploaded '+ date +' by '+ user;
        break;
      default:

        //none
        break;
    }

    return (

      <span className='light-text'>

        { status }
      </span>
    );
  }

  render(){
    
    var fileUpload = this['props']['fileUpload'];

    return (

      <div className='list-entry-group'>

        <div className='list-entry-icon'>

          <span className='glyphicon glyphicon-file'></span>
        </div>
        <div className='list-entry-title-group'>
          
          <div className='list-entry-file-name'>
            
            { this.parseFileName() }
          </div>
          <div className='list-entry-status'>

            { this.parseFileStatus() }
          </div>
        </div>
        <UploadMeter fileUpload={ fileUpload } />
        <UploadControl fileUpload={ fileUpload } />
      </div>
    );
  }
}