import * as React from 'react';

import UploadControl from './upload-control';

export default class ListUploadFailed extends React.Component{

  constructor(props){

    super(props);
  }

  removeFile(){

    this['props']['callbacks'].removeFailed(this['props']['failedFile']);
  }

  selectUpload(){

    /*
     * We are using a button to surragate the file input so we may have 
     * a custom browse button.
     */
    var failedFile = this['props']['failedFile'];
    document.getElementById('file-selector-'+ failedFile['reactKey']).click();
  }

  fileSelected(event){

    if(event === undefined) return;
    var fileSelector = event.target;
    var file = fileSelector.files[0];
    this['props']['callbacks'].recoverUpload(file, this['props']['failedFile']);
  }

  render(){

    var failedFile = this['props']['failedFile'];
    var currentFileSize = PARSE_BYTES(failedFile['currentBytePosition'], 1000);

    var uploadControl = {

      'fileUpload': this['props']['failedFile'],
      'callbacks': {

        'selectUpload': this['selectUpload'].bind(this),
        'cancelUpload': this['removeFile'].bind(this)
      }
    };

    var fileSelector = {

      'id': 'file-selector-'+ failedFile['reactKey'],
      'className': 'file-selector',
      'type': 'file',
      'name': 'upload-file-'+ failedFile['reactKey'],
      'onChange': this['fileSelected'].bind(this)
    };

    return (

      <tr className='list-entry-group'>

        <td className='list-entry-icon'>
        </td>
        <td className='list-entry-title-group'>

          <div className='list-entry-file-name' title={failedFile['fileName']}>

            {failedFile['fileName']}
          </div>
          <div className='list-entry-status' title='The current file status.'>

            <span className='light-text'>

              {'attempted at: '+ PARSE_TIMESTAMP(failedFile['startTimestamp'])}
            </span>
          </div>
        </td>
        <td className='list-entry-project-group'>

          <div className='list-entry-project-name'>

            { failedFile['projectName'] }
          </div>
          <div className='list-entry-role'>

            <span className='light-text'>

              { failedFile['role'] }
            </span>
          </div>
        </td>
        <td className='list-entry-title-group'>

          <div className='dark-text' style={{ 'fontWeight': 900 }} >

            { 'current data uploaded: '+ currentFileSize }
          </div>
          <div className='list-entry-role'>
            
            <span className='light-text'>

              { 'this file upload was interrupted.' }
            </span>
          </div>
        </td>
        <input { ...fileSelector } />
        <UploadControl { ...uploadControl } />
      </tr>
    )
  }
}