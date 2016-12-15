import * as ReactRedux from 'react-redux';
import AdminUI from './admin-ui';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    appState: state['appState']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {};
}

const AdminUIContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(AdminUI);

export default AdminUIContainer;