import * as ReactRedux from 'react-redux';
import UploadControl from './upload-control';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {};
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    startUpload: (file)=>{

      var action = { type: 'START_UPLOAD', data: file };
      dispatch(action);
    },

    pauseUpload: (file)=>{

      var action = { type: 'PAUSE_UPLOAD', data: file };
      dispatch(action);
    },

    cancelUpload: (file)=>{

      var action = { type: 'CANCEL_UPLOAD', data: file };
      dispatch(action);
    }
  };
}

const UploadControlContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(UploadControl);

export default UploadControlContainer;