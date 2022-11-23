import React, { useCallback } from 'react';
import {connect} from 'react-redux';

function UploadDialog({ dismissDialog, startFileUpload, startDirectoryUpload }) {
  const onUploadDirectory = useCallback(() => {
    dismissDialog();
    startDirectoryUpload();
  }, [startDirectoryUpload, dismissDialog]);

  const onUploadFile = useCallback(() => {
    dismissDialog();
    startFileUpload();
  }, [startFileUpload, dismissDialog]);

  return (<div className='upload-dialog'>
    <div className='title'>Which type of upload?</div>
    <div className=''>
      <button className='button' onClick={onUploadDirectory}>
        Upload Directory
      </button>
      <button className='button' onClick={onUploadFile}>
        Upload File
      </button>
    </div>
  </div>);
}

export default connect(
  null,
  (dispatch) => ({
    dismissDialog: () => dispatch({type:'DISMISS_DIALOG'})
  })
)(UploadDialog);
