import * as ReactRedux from 'react-redux';
import PermissionEdit from './permission-edit';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    adminInfo: state['adminInfo']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    'addPermission': ()=>{

      var action = { 'type': 'ADD_PERMISSION' };
      dispatch(action);
    },

    'savePermission': (permission)=>{

      var action = { 'type': 'SAVE_PERMISSION', 'permission': permission };
      dispatch(action);
    },

    'removeUnsavedPermission': (reactKey)=>{

      var action = { 

        'type': 'REMOVE_UNSAVED_PERMISSION', 
        'reactKey': reactKey 
      };
      dispatch(action);
    },

    'downloadPermissions': ()=>{

      var action = { 'type': 'DOWNLOAD_PERMISSIONS' };
      dispatch(action);
    },

    'uploadPermissions': (file)=>{

      var action = { 'type': 'UPLOAD_PERMISSIONS', 'file': file };
      dispatch(action);
    }
  };
}

const PermissionEditContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(PermissionEdit);

export default PermissionEditContainer;