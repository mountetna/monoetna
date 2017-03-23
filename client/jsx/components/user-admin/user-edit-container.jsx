import * as ReactRedux from 'react-redux';
import UserEdit from './user-edit';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    'adminInfo': state['adminInfo']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {

    'logoutAll': ()=>{

      var action = { 'type': 'LOGOUT_ALL' };
      dispatch(action);
    }
  }
}

const UserEditContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(UserEdit);

export default UserEditContainer;