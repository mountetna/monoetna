import * as React from 'react'

import UploadMeter from './upload-meter'
import UploadControl from './upload-control'

export default class ListUpload extends React.Component{

  constructor(){

    super();
  }

  parseFileName(){

    var file = this['props']['fileUpload'];
    var origName = file['originalName'];
    var fileName = file['fileName'];

    /*
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
    */

    return fileName;
  }

  parseFileStatus(){

    var file = this['props']['fileUpload'];
    var status = file['status'];
    var date = PARSE_TIMESTAMP(Date.now() / 1000);
    var user = file['userEmail'];
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

      <tr className='list-entry-group'>

        <td className='list-entry-icon'>
        </td>
        <td className='list-entry-title-group'>

          <div className='list-entry-file-name'>
            
            { this.parseFileName() }
          </div>
          <div className='list-entry-status'>

            { this.parseFileStatus() }
          </div>
        </td>
        <td className='list-entry-project-group'>

          {/*
          <div className='list-entry-file-name' title='The projec this file belongs to.'>

            { 'prjkt' }
          </div>
          <div className='list-entry-status' title='Your project permission for this file.'>

            <span className='light-text'>
              
              { 'administrator' }
            </span>
          </div>
        */}
        </td>
        <UploadMeter fileUpload={ fileUpload } />
        <UploadControl fileUpload={ fileUpload } />
      </tr>
    );
  }
}