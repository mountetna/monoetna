import * as ReactRedux from 'react-redux';
import PermissionEdit from './permission-edit';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    appState: state['appState']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    addPermission: ()=>{

      var action = { 'type': 'ADD_PERMISSION' };
      dispatch(action);
    },

    savePermission: (permission)=>{

      var action = { 'type': 'SAVE_PERMISSION', 'permission': permission };
      dispatch(action);
    },

    removeUnsavedPermission: (permissionId)=>{

      var action = { 

        'type': 'REMOVE_UNSAVED_PERMISSION', 
        'permissionId': permissionId 
      };
      dispatch(action);
    }
  };
}

const PermissionEditContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(PermissionEdit);

export default PermissionEditContainer;