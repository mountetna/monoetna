import * as ReactRedux from 'react-redux';
import ListBody from './list-body';

const mapStateToProps = (state, ownProps)=>{

  // Blend the specific user data permissions in with the file data.
  var permissions = state['userInfo']['permissions'];

  var permissionMapper = (listElem, perms)=>{

    for(var c = 0; c < perms['length']; ++c){

      if(parseInt(listElem['projectId']) == perms[c]['projectId']){

        listElem['projectName'] = perms[c]['projectName'];
        listElem['role'] = perms[c]['role'];
        listElem['groupId'] = perms[c]['groupId'];
        break;
      }
    }
  }

  for(var a = 0; a < state['fileData']['fileList']['length']; ++a){

    permissionMapper(state['fileData']['fileList'][a], permissions);
  }

  for(var b = 0; b < state['fileData']['fileFails']['length']; ++b){

    permissionMapper(state['fileData']['fileFails'][b], permissions);
    if(state['fileData']['fileFails'][b]['status'] != 'cancelled'){

      state['fileData']['fileFails'][b]['status'] = 'failed';
    }
  }

  // state == redux store
  return {

    'userInfo': state['userInfo'],
    'fileData': state['fileData']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    initializeUpload: (file)=>{

      var action = { 'type': 'AUTHORIZE_FILE', 'uploadFile': file };
      dispatch(action);
    },

    queueUpload: (redisIndex)=>{

      var action = { 'type': 'QUEUE_UPLOAD', 'redisIndex': redisIndex };
      dispatch(action);
    },

    pauseUpload: (redisIndex)=>{

      var action = { 'type': 'PAUSE_UPLOAD', 'redisIndex': redisIndex };
      dispatch(action);
    },

    cancelUpload: (redisIndex)=>{

      var action = { 'type': 'CANCEL_UPLOAD', 'redisIndex': redisIndex };
      dispatch(action);
    },

    recoverUpload: (file, fileMetadata)=>{

      var action = { 

        'type': 'RECOVER_UPLOAD',
        'uploadFile': file,
        'fileMetadata': fileMetadata
      };

      dispatch(action);
    },

    removeFile: (fileMetadata)=>{

      var action = { 'type': 'REMOVE_FILE', 'fileMetadata': fileMetadata };
      dispatch(action);
    }
  };
}

const ListBodyContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(ListBody);

export default ListBodyContainer;