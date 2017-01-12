import * as ReactRedux from 'react-redux';
import ListBody from './list-body';

const mapStateToProps = (state, ownProps)=>{

  // Blend the specific user data permissions in with the file data.
  var fileList = state['fileData']['fileList'];
  var permissions = state['userInfo']['permissions'];

  for(var a = 0; a < fileList['length']; ++a){

    for(var b = 0; b < permissions['length']; ++b){

      if(parseInt(fileList[a]['projectId']) == permissions[b]['projectId']){

        fileList[a]['projectName'] = permissions[b]['projectName'];
        fileList[a]['role'] = permissions[b]['role'];
        fileList[a]['groupId'] = permissions[b]['groupId'];
        break;
      }
    }
  }

  state['fileData']['fileList'] = fileList;

  // state == redux store
  return {

    'userInfo': state['userInfo'],
    'fileData': state['fileData']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    authorizeFile: (file)=>{

      var action = { type: 'AUTHORIZE_FILE', data: file };
      dispatch(action);
    }
  };
}

const ListBodyContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(ListBody);

export default ListBodyContainer;