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

      return { 'width': '0%' };
    }
    else{

      return { 'width': String((bytesUploaded/fileSize)*100) + '%' };
    }
  }

  parseUploadBytes(){

    var file = this['props']['fileUpload'];
    var fileSize = PARSE_BYTES(file['fileSize'], true);
    var bytesUploaded = PARSE_BYTES(file['currentBytePosition'], true);

    var uploadInfoProps = {

      'className': 'upload-meter-info light-text',
      'style': { 'float': 'left' }
    }

    var bytesUploadedProps = {

      'className': 'dark-text',
      'style': { fontWeight: 900 },
      'title': 'kiloBYTES uploaded.'
    }

    return (

      <div { ...uploadInfoProps }>

          <span { ...bytesUploadedProps }>

            { bytesUploaded }
          </span>
          { ' of '+ fileSize +' uploaded'}
      </div>
    );
  }

  parseUploadSpeed(){

    var file = this['props']['fileUpload'];
    if('uploadSpeed' in file){

      if(isNaN(file['uploadSpeed'])) return '';

      var speed = (file['uploadSpeed']/1000).toFixed(1);

      var bitSpeedProps = {

        'className': 'dark-text',
        'style': { fontWeight: 900 },
        'title': 'kiloBITS per second.'
      }

      return (

        <div className='upload-meter-info' style={{ 'float': 'right' }}>

          <span { ...bitSpeedProps }>
            
            { speed + ' kbps' }
          </span>
        </div>
      );
    }
    else{

      return '';
    }
  }

  render(){
    
    return (

      <td className='upload-meter-group'>
 
        <div className='upload-meter-tray'>

          <div className='upload-meter-bar' style={ this.calcUploadPercent() }>
          </div>
        </div>

        { this.parseUploadBytes() }
        { this.parseUploadSpeed() }
      </td>
    );
  }
}