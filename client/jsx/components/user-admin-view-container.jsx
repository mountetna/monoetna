import * as ReactRedux from 'react-redux';
import UserAdminView from './user-admin-view';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    'userInfo': state['userInfo']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {};
}

const UserAdminViewContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(UserAdminView);

export default UserAdminViewContainer;