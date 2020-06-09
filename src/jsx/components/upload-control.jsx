import * as React from 'react';
import { connect } from 'react-redux';

const UploadButton = ({onClick, icon}) =>
  <button className='upload-control-btn' onClick={onClick}>
    <span className={`fas fa-fw fa-${icon}`}></span>
  </button>;

const UploadControl = ({ upload, continueUpload, pauseUpload, cancelUpload, selectUpload }) => {
  let invoke = (callback) => () => callback(upload);

  let buttonProps = {
    paused: { icon: 'play', onClick: invoke(continueUpload) },
    active: { icon: 'pause', onClick: invoke(pauseUpload) },
    failed: { icon: 'retweet', onClick: invoke(selectUpload) }
  };

  return (
    <div className='upload-control-group'>
      <UploadButton { ...buttonProps[upload.status] }/>
      <UploadButton icon='times' onClick={invoke(cancelUpload)}/>
    </div>
  );
};

export default connect(
  // map state
  null,

  // map dispatch
  null
)(UploadControl);
