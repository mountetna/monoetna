import * as ReactRedux from 'react-redux';
import UserEdit from './user-edit';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    'adminInfo': state['adminInfo']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {};
}

const UserEditContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(UserEdit);

export default UserEditContainer;