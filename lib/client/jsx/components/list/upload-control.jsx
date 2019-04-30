import * as React from 'react';
import { connect } from 'react-redux';

const UploadButton = ({onClick, icon}) =>
  <button className='upload-control-btn' onClick={onClick}>
    <span className={`fas fa-fw fa-${icon}`}></span>
  </button>;

class UploadControl extends React.Component{
  render(){
    let { upload, queueUpload, pauseUpload, cancelUpload, selectUpload } = this.props;
    let invoke = (callback) => () => callback(upload);

    let buttonProps = {
      paused: { icon: 'play', onClick: invoke(queueUpload) },
      active: { icon: 'pause', onClick: invoke(pauseUpload) },
      failed: { icon: 'retweet', onClick: invoke(selectUpload) }
    };

    return (
      <div className='upload-control-group'>
        <UploadButton { ...buttonProps[upload.status] }/>
        <UploadButton icon='times' onClick={invoke(cancelUpload)}/>
      </div>
    );
  }
}

export default connect(
  // map state
  null,

  // map dispatch
  (dispatch) => ({
    cancelUpload: (upload) => dispatch({ type: 'CANCEL_UPLOAD', upload }),
    queueUpload: (upload) => dispatch({ type: 'QUEUE_UPLOAD', upload }),
    pauseUpload: (upload) => dispatch({ type: 'PAUSE_UPLOAD', upload }),
  })
)(UploadControl);
