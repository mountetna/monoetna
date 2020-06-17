import * as React from 'react';
import { byteFormat } from '../utils/format';
import Icon from './icon';

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

const ICONS={
  complete: 'check',
  queued: 'clock',
  active: 'upload',
  paused: 'pause'
};

const UploadStatus = ({upload: {status}}) =>
  <div className='upload-meter-status'>
    { ICONS[status] ? <Icon icon={ICONS[status]} className={status}/> : status }
  </div>;

const UploadCompletion = ({upload}) => {
  let { file_size, current_byte_position } = upload;

  return <div className='upload-meter-completion'>
    <span className='completed' title='kilobytes uploaded'>
      { byteFormat(current_byte_position, true) }
    </span> of { byteFormat(upload.file_size, true) }
  </div>
}

const UploadSpeed = ({upload}) => {
  let { paused, upload_speeds } = upload;

  if (!upload_speeds.length || paused) return null;

  let avgSpeed = upload_speeds.reduce((sum,t)=>sum+t, 0) / upload_speeds.length;

  // rough calculation of the upload speed in kbps for the UI
  let upload_speed = avgSpeed * 8 * 1000;

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
        <UploadStatus upload={ upload }/>
        <UploadCompletion upload={ upload }/>
        <UploadSpeed upload={ upload }/>
      </div>
    );
  }
}
