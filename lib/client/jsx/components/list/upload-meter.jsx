import * as React from 'react';
import { byteFormat } from '../../utils/format';

const UploadBar = ({upload}) => {
  let { file_size, current_byte_position } = upload;
  return <div className='upload-meter-tray'>
    <div className='upload-meter-bar' style={
      {
        width: file_size == 0 ? '100%' : `${(current_byte_position/file_size)*100}%`
      }
    }/>
  </div>
}

const UploadCompletion = ({upload}) => {
  let { file_size, current_byte_position } = upload;

  return <div className='upload-meter-completion'>
    <span className='completed' title='kilobytes uploaded'>
      { byteFormat(current_byte_position, true) }
    </span> of { byteFormat(upload.file_size, true) }
    { file_size == current_byte_position && <i className='done fas fa-check'/> }
  </div>
}

const UploadSpeed = ({upload}) => {
  let { paused, upload_speed } = upload;

  if (!upload_speed || isNaN(upload_speed) || paused) return null;

  return (
    <div className='upload-meter-speed'>
      <span title='kilobits per second'>
        { byteFormat(upload_speed, false, 'bps') }
      </span>
    </div>
  );
}

export default class UploadMeter extends React.Component{
  render(){
    let { upload } = this.props;
    return (
      <div className='upload-meter-group'>
        <UploadBar upload={ upload }/>
        <UploadCompletion upload={ upload }/>
        <UploadSpeed upload={ upload }/>
      </div>
    );
  }
}
