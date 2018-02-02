import * as ReactRedux from 'react-redux';
import MetisUI from './metis-ui';

const mapStateToProps = (state, ownProps)=>{
  // state == redux store
  return {
    userInfo: state.userInfo,
    fileData: state.fileData
  };
}

const mapDispatchToProps = (dispatch, ownProps)=>{
  return {};
}

const MetisUIContainer = ReactRedux.connect(
  mapStateToProps,
  mapDispatchToProps,
)(MetisUI);

export default MetisUIContainer;
