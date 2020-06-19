import * as React from 'react';
import { connect } from 'react-redux';

import {
  CANCEL_UPLOAD,
  CONTINUE_UPLOAD,
  PAUSE_UPLOAD
} from '../actions/upload_actions';

const UploadButton = ({ onClick, icon }) => (
  <button className="upload-control-btn" onClick={onClick}>
    <span className={`fas fa-fw fa-${icon}`}></span>
  </button>
);

const UploadControl = ({ upload, continueUpload, pauseUpload, cancelUpload, selectUpload }) => {
  let invoke = (callback) => () => callback(upload);

  let buttonProps = {
    paused: { icon: 'play', onClick: invoke(continueUpload) },
    active: { icon: 'pause', onClick: invoke(pauseUpload) },
    failed: { icon: 'retweet', onClick: invoke(selectUpload) }
  };

  return (
    <div className='upload-control-group'>
      <UploadButton { ...buttonProps[upload.status] } />
      <UploadButton icon='times' onClick={invoke(cancelUpload)} />
    </div>
  );
};

export default connect(
  // map state
  null,

  // map dispatch
  (dispatch) => ({
    cancelUpload: (upload) => dispatch({ type: CANCEL_UPLOAD, upload }),
    continueUpload: (upload) => dispatch({ type: CONTINUE_UPLOAD, upload }),
    pauseUpload: (upload) => dispatch({ type: PAUSE_UPLOAD, upload })
  })
)(UploadControl);
