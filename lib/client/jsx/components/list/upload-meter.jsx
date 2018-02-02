import * as React from 'react';
import { byteFormat, dateFormat } from '../../utils/format';

export default class UploadMeter extends React.Component{
  constructor(){
    super();
  }

  calcUploadPercent(){
    let { fileUpload } = this.props;
    let { fileSize, currentBytePosition } = fileUpload;

    if (fileSize == 0) {
      return { width: '0%' };
    }
    else {
      return { width: String((currentBytePosition/fileSize)*100) + '%' };
    }
  }

  parseUploadBytes(){
    let { fileUpload } = this.props;
    let { fileSize, currentBytePosition } = fileUpload;

    fileSize = byteFormat(fileSize, true);
    let bytesUploaded = byteFormat(currentBytePosition, true);

    let uploadInfoProps = {
      className: 'upload-meter-info light-text',
      style: { float: 'left' }
    }

    let bytesUploadedProps = {
      className: 'dark-text',
      style: { fontWeight: 900 },
      title: 'kiloBYTES uploaded.'
    }

    return (
      <div { ...uploadInfoProps }>
          <span { ...bytesUploadedProps }>
            { bytesUploaded }
          </span>
          { ` of ${fileSize} uploaded`}
      </div>
    );
  }

  parseUploadSpeed(){
    let file = this.props.fileUpload;

    if('uploadSpeed' in file && file.status == 'active'){
      if(isNaN(file.uploadSpeed)) return '';

      let speed = byteFormat(file.uploadSpeed, 1024, true);

      let bitSpeedProps = {
        className: 'dark-text',
        style: { fontWeight: 900 },
        title: 'kiloBITS per second.'
      }

      return (
        <div className='upload-meter-info' style={{ float: 'right' }}>
          <span { ...bitSpeedProps }>
            { speed }
          </span>
        </div>
      );
    }
    else {
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
