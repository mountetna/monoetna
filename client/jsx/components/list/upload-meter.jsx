import * as React from 'react'

export default class UploadMeater extends React.Component{

  constructor(){

    super();
  }

  calcUploadPercent(){

    var file = this['props']['fileUpload'];
    var fileSize = file['fileSize'];
    var bytesUploaded = file['currentBytePosition'];

    if(fileSize == 0){

      return { width: '0%' };
    }
    else{

      return { width: String((bytesUploaded/fileSize)*100) + '%' };
    }
  }

  parseUploadBytes(){

    var file = this['props']['fileUpload'];
    var fileSize = PARSE_BYTES(file['fileSize'], true);
    var bytesUploaded = PARSE_BYTES(file['currentBytePosition'], true);

    return (

      <div className='list-entry-status light-text'> 

        { bytesUploaded +" of "+ fileSize }
      </div>
    );
  }

  render(){
    
    return (

      <div className='upload-meter-group'>

        <div className='upload-meter-tray'>

          <div className='upload-meter-bar' style={ this.calcUploadPercent() }>
          </div>
        </div>
        <div className='upload-meter-info'>

          { this.parseUploadBytes() }
        </div>
      </div>
    );
  }
}