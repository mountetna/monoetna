import * as ReactRedux from 'react-redux';
import ListBody from './list-body';

const mapStateToProps = (state, ownProps)=>{

  // state == redux store
  return {

    fileList: state['appState']['fileList'],
    fileUploads: state['appState']['fileUploads'],
    permissions: state['appState']['userInfo']['permissions']
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{

  return {};
}

const ListBodyContainer = ReactRedux.connect(

  mapStateToProps,
  mapDispatchToProps,
)(ListBody);

export default ListBodyContainer;